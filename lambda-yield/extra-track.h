// Extra-track helpers shared by ntrack-class.cc and yield-compare.cc.
#ifndef LAMBDA_YIELD_EXTRA_TRACK_H
#define LAMBDA_YIELD_EXTRA_TRACK_H
#include <algorithm>
#include <cmath>
#include <vector>

namespace extra_track {
constexpr double kZTarget = -143.;
// A good extra track has more than kShortHits hits and is not an upstream beam-window track.
constexpr size_t kShortHits = 6;
// A good extra track closer than kVertexTrackDca to the production vertex vetoes the event
// in the vertex ntrack: it comes from the reaction, the farther ones are accidental.
constexpr double kVertexTrackDca = 30.; // mm

// Same test as IsUpstreamBeamWindowTrack of DstTPCLambdaEta.cc.
inline bool IsUpstream(const std::vector<double> &x, const std::vector<double> &z)
{
  bool any = false;
  for (size_t i = 0; i < std::min(x.size(), z.size()); ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(z[i])) continue;
    any = true;
    if (z[i] >= kZTarget || x[i] <= -50. || x[i] >= 50.) return false;
  }
  return any;
}
// Closest distance of a TPC helix to a global position. The helix is in the TPC
// local frame, (x, y, z)_local = (-x, z - Z_T, y)_global, as LocalToGlobal of
// TPCLocalTrackHelix.cc; t is scanned over the hit range with a margin of pi/2.
inline double HelixDistance(double cx, double cy, double z0, double r, double dz, double tmin, double tmax,
                           double gx, double gy, double gz)
{
  const double x = -gx, y = gz - kZTarget, z = gy;
  auto d2 = [&](double t) {
    const double ex = cx + r * cos(t) - x, ey = cy + r * sin(t) - y, ez = z0 + dz * r * t - z;
    return ex * ex + ey * ey + ez * ez;
  };
  const double lo = std::min(tmin, tmax) - 0.5 * M_PI, hi = std::max(tmin, tmax) + 0.5 * M_PI;
  const int n = 400;
  int best = 0;
  for (int i = 1; i <= n; ++i)
    if (d2(lo + (hi - lo) * i / n) < d2(lo + (hi - lo) * best / n)) best = i;
  double a = lo + (hi - lo) * std::max(best - 1, 0) / n, b = lo + (hi - lo) * std::min(best + 1, n) / n;
  const double g = 0.5 * (sqrt(5.) - 1);
  for (int it = 0; it < 60; ++it) {
    const double c = b - g * (b - a), e = a + g * (b - a);
    if (d2(c) < d2(e)) b = e;
    else a = c;
  }
  return sqrt(d2(0.5 * (a + b)));
}
} // namespace extra_track
#endif
