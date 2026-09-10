# Sensitivity 매크로 계보 및 이관 자료

조사일: 2026-09-10. 원본: `/Users/ihaein/Work/E72/Simul/tpc_simple_analyzer/macro`.
아래 경로는 별도 표시가 없으면 원본 macro 기준이다. 원본 코드·데이터는 수정하거나 이동하지 않았다.
코드 비교와 텍스트 데이터 검사를 수행했으며 ROOT 매크로 실행 및 ROOT 내부 객체 검증은 수행하지 않았다.

## 1. 버전 순서

사용자 설명에 따른 작업 계보는 **polarity 초기 작업 → macro로 복사하여 빔타임 중 갱신**이다.
다만 현재 polarity도 후속 수정되어 초기 스냅샷 그대로는 아니다. 원본 디렉터리는 Git 저장소가 아니므로 정확한 복사 날짜와 모든 중간 버전을 복원할 수 없다. 아래 날짜는 파일 수정 시각이며 작성/실행 시각의 증거는 아니다.

| 단계/시점 | 파일 | 현재 남아 있는 내용과 의미 |
|---|---|---|
| 기초 도구, 2025-03-24~31 | `diff_cross.cc`, `acceptance.cc`, `total_E72.cc` | 상위 파일의 수정 시각은 polarity보다 오래됨. MC 각도분포·acceptance와 문헌 total cross section 도구가 이미 존재했음을 보여줌 |
| 직접 편극 계산, 06-25 | `polar.cc` | CSFlat MC의 입자 4-vector와 Λ 붕괴 양성자 방향에서 편극 계산. 뒤의 통계 민감도 예측과 별도 도구 |
| polarity 작업, 09-28~10-08 | `polarity/{diff_cross,acceptance,polar_th}.cc`, `README` | 724~770 MeV/c의 연속 2 MeV/c 간격 24-bin 체계. README에 계산 순서가 남아 있음 |
| 빔 프로파일 기반 계획, 10-09 | `polarity/make_beamfile.cc` | run344 프로파일로 총 7일 중 고정빔 3.5~5.5일의 다섯 시나리오를 생성 |
| 예상 빔 기반 yield, 11-20 | `yield_E72.cc` | 가우시안 빔을 직접 생성. trigger efficiency 0.78. run69 ROOT 파일을 열지만 실제 yield 계산에는 사용하지 않음 |
| 빔 수 입력 교체, 11-23 | `data/convert_csv.cc`, `yield.cc`, `polarity/yield.cc` | 두 yield 모두 `Nbeam_251123.root:hMom2`를 사용하고 trigger efficiency 0.5. 상위는 15-bin, polarity는 24-bin. 수정 시각은 각각 08:14와 08:17이므로 polarity에도 후속 반영됨 |
| 단면적 예측 수정, 12-06~07 | 두 `crystalball.cc` | 상위는 15-bin, polarity는 24-bin으로 선언. 일부 fit/그림 코드 활성화 상태가 다름. polarity 쪽은 12-07 수정되어 초기 원본으로 볼 수 없음 |
| 상위 편극 예측 수정, 12-07 | `polar_th.cc` | 15-bin. 735 MeV/c 모델 비교의 오차를 첫 10개 운동량 bin 합산 통계로 변경 |

따라서 두 폴더는 **24-bin 계열과 15-bin 계열**로 보존하고, `yield_E72.cc → yield.cc`를 빔 입력 방식 변경으로 따로 기록하는 편이 정확하다. 단순히 폴더 전체를 v1/v2로 부르면 실제 수정 이력이 가려진다.

## 2. 계산 순서와 역할

```text
Nbeam_251123.csv → data/convert_csv.cc → Nbeam_251123.root (hMom2)
                                             ↓
                                          yield.cc
                              ┌──────────────┴─────────────┐
                        entry_mom.txt              luminosity_mom.txt
                              ↓                            │
CSOn MC ───────────────→ diff_cross.cc                      │
                              ↓                            │
                       entry_mom_cos.txt ───────┬───────────┤
                              ↓                │           │
                         polar_th.cc           ↓           ↓
Acceptance_study MC → acceptance.cc ───────→ crystalball.cc
                         acceptance_mom_cos.txt
```

실행 순서는 각 버전 폴더를 작업 디렉터리로 하여 `yield.cc → diff_cross.cc → polar_th.cc`; `crystalball.cc`에는 추가로 `acceptance.cc` 출력이 필요하다. ROOT 상대경로는 매크로 파일 위치가 아니라 실행 작업 디렉터리에 의존한다. 현 코드에는 아래 오류 후보가 있어 이 순서는 의존관계 설명이며 실행 검증 완료를 뜻하지 않는다.

| 매크로 | 입력 및 선택 조건 | 산출물/계산 |
|---|---|---|
| `yield.cc` 양쪽 | 빔 수 `hMom2`, 코드 내 total cross section 배열 | `param/entry_mom.txt`, `param/luminosity_mom.txt` |
| `yield_E72.cc` | 가우시안 예상 빔; 685/705/725/745/765 각각 0.5일 + 735에서 5.5일 | 같은 이름의 두 param 파일을 출력하므로 현재 yield 결과를 덮어쓸 수 있음 |
| `diff_cross.cc` 양쪽 | CSOn `g4hyptpc_light`; `mom_kaon_lab`, `cos_theta`, `trig_flag`, `decay_particle_code`; 2112 제외, trigger nonzero | 통과 MC의 각도별 비율에 entry_mom을 배분 → `param/entry_mom_cos.txt`. 이 파일 자체는 differential cross section 값이 아니라 예상 사건 수 |
| `acceptance.cc` 양쪽 | Acceptance_study의 같은 tree/branches; 2112 제외 후 trigger nonzero 비율 | `param/acceptance_mom_cos.txt` |
| `polar_th.cc` 양쪽 | entry_mom_cos + 하드코딩된 Crystal Ball 편극 두 세트 및 D03/P03 모델 좌표 | spline 모델 곡선과 통계오차 `sqrt(3)/(0.65*sqrt(N))`; 운동량별 통계오차 |
| `crystalball.cc` 양쪽 | luminosity, entry_mom_cos, acceptance_mom_cos + 하드코딩된 15 운동량 × 9 각도 단면적/오차 | 단면적 fit 및 예상 오차 `sqrt(N)/(luminosity * Δcos * 2π * acceptance)` |
| `polar.cc` | CSFlat `g4hyptpc`; BEAM/PRM/SEC, cos_theta; Λ→pπ 붕괴 입자 | 9 각도 bin, 720~770의 10 운동량 bin; 방향 평균 기반 편극 |
| `total_E72.cc` | 코드 내 문헌 total cross section 배열 | 문헌 단면적 비교 그림 |

상위 15-bin 중심: `724,726,728,730,732,734,738,742,746,750,754,758,762,766,770`.
polarity 24-bin 중심: `724,726,...,770`. 둘 다 각 중심 ±1 MeV/c, 각도는 15등분이다. 상위 15-bin은 734 이후 운동량 구간에 빈틈이 있다.

현재 param도 상위는 15/225행, polarity는 24/360행으로 확인했다. 일부 param의 수정일은 2026-09-10이므로 2025 빔타임 당시 결과의 보존본이라고 단정할 수 없다.

## 3. 사용하는 데이터와 출처

### 현재 sensitivity 계산의 필수 입력

| 실제 위치 (analyzer 기준) | 용도 | 파일 존재/크기 |
|---|---|---|
| `data/7201_7202_CSOn.root` | 양쪽 diff_cross의 MC 각도 shape | 확인, 75,131,532 bytes |
| `data/7201_7202_Acceptance_study.root` | 양쪽 acceptance의 MC 효율 | 확인, 75,100,937 bytes |
| `macro/data/Nbeam_251123.root` | 양쪽 yield의 빔 수 `hMom2` | 확인, 9,041 bytes |
| `macro/data/Nbeam_251123.csv` | 바로 위 ROOT의 변환 원본 | 확인, 856행; 중심 575.25~1002.75, 간격 0.5 MeV/c |
| 두 버전의 `param/*.txt` | 중간 결과 및 재현용 현 상태 | 각 버전 네 파일 모두 확인 |

`convert_csv.cc`는 CSV 2열을 bin content로 넣고 4개씩 합쳐 2 MeV/c 히스토그램 `hMom2`를 만든다. 코드가 의도하는 범위는 575~1003, 214 bins이며 yield 설정과 일치한다. CSV가 어떤 run 목록, scaler 보정, beam cut, DAQ 보정을 거쳐 만들어졌는지는 이 폴더의 코드만으로 확인되지 않는다. 파일명 날짜를 실제 수집기간으로 해석해서는 안 된다.

문헌 total cross section은 ROOT 입력이 아니라 코드 안 배열이다. 주석의 DOI는 `10.1007/BF02785525` (27점, 777~1842)와 `10.1103/PhysRevC.64.055205` (17점, 720~770). yield는 이들을 합치고 700/718에 0을 추가한 뒤 `TGraph::Eval`로 보간한다. 이는 코드에 기록된 출처이며 논문 원문 대조는 하지 않았다. crystalball의 differential 자료와 polar_th의 편극 자료는 Crystal Ball로 표시되어 있으나 모델 좌표를 추출한 원문·digitization 기록은 확인되지 않았다.

### 이전·보조 데이터

| 파일 (`macro/data` 기준, 예외 명시) | 현재 사용 상태 |
|---|---|
| `beam_profile_run344_layer13_center.root` | make_beamfile의 `tr` tree, pInx/pIny/pInz를 사용; 크기에서 운동량을 구해 ×1000 |
| `n_kaon_scan.root` | make_beamfile 출력. hist_beam_35/40/45/50/55. 현재 두 yield의 입력은 아님 |
| `n_kaon.root` | yield_E72의 mom_all 로딩 코드는 주석 처리됨 |
| `beam.k.run69_0130.root` | yield_E72가 열기만 하고 histogram/tree를 읽지 않음. yield의 실제 beam 입력으로 분류하면 안 됨 |
| `run00344_Hodoscope.root`, `run00344_BcOutTracking_layer12.root` | 조사한 sensitivity 매크로에서 직접 읽지 않음. 프로파일 생성과의 연결은 추가 자료 필요 |
| `e72_defocus_HS_off.root` | e72_defocus.cc의 TPC hit 그림; 별도 geometry/defocus 확인용 |
| `e72_defocus_e72_defocus_check.root` | 조사한 sensitivity 매크로의 활성 입력 아님 |
| analyzer `data/7201_7202_CSFlat.root` | polar.cc 전용, 존재 확인, 2,699,573,882 bytes. 주 sensitivity 파이프라인에는 불필요 |

make_beamfile은 총 7일, 고정빔 3.5~5.5일/나머지 scan, beam power 90/82, target beam fraction 0.69, BAC-target 0.35 m의 붕괴 보정, 735 MeV/c intensity를 사용한다. scan은 620~950 범위에 균등하게 배분한다. 이 시나리오가 Nbeam CSV 생성에 사용됐다는 연결은 코드에서 확인되지 않는다.

## 4. 이관 전 확인할 계산상 문제

1. **상위 15-bin의 운동량 빈틈**: diff_cross는 전체 723~771 범위만 먼저 자른다. 이후 어느 ±1 구간에도 속하지 않으면 mombin=-9999인 채 배열을 접근한다. acceptance도 mombin 미초기화 상태가 가능하다. 두 계열 모두 최상단 771 경계 처리도 확인 필요.
2. **polarity/crystalball의 24/15 불일치**: n_mombin=24인데 poly 및 문헌 운동량 배열은 15개이고 g_diff_exp도 15개로 선언되어 있다. 24회 계산 루프에서 범위를 넘는다. 단순히 최신 파일을 실행해서 비교할 수 있는 상태가 아니다.
3. **편극 운동량 label/통계의 불일치**: polarity의 [5]/[13]은 734/750인데 모델 label은 735/765다. 상위 [13]은 766이다. 상위 첫 그림의 통계는 첫 10개 중심(724~750)을 합친 값이므로 735 단일 bin의 민감도로 해석하면 안 된다.
4. **정규화 확인**: yield의 출력 entry는 양 끝의 yield 평균을 구간폭으로 다시 나눈다. 입력이 2 MeV/c당 빔 수인 상태에서 후속 코드가 이것을 사건 수로 쓰므로 bin당 수/밀도 정의를 맞춰야 한다. luminosity는 E_eff만 포함하고 yield는 trigger_eff와 Λ 분기비도 포함한다. crystalball 오차와 함께 효율 정의를 검토해야 한다.
5. **fit 함수**: 두 crystalball의 Legendre2에서 `3/2`, `1/2`는 C++ 정수 나눗셈이어서 의도한 P2와 다르다.
6. **acceptance 오차**: acceptance가 이미 0~1 비율인데 오차 계산에서 다시 0.01을 곱한다. 저장 TXT에는 오차가 없지만 그림의 오차에는 영향을 준다.
7. **각도 표기**: diff_cross는 cos_L=-cos_theta를 만들고 사용하지 않으면서 축에는 Λ를 표시한다. acceptance/polar_th는 η 표기다. polar_th의 문헌 마지막 cos 중심도 0.85로 적혀 있어 원자료 대조가 필요하다.

현재 yield의 공통 수치는 target density 0.07085 g/cm³, 길이 7.139 cm, 원자량 1.008, E_eff=0.8, trigger_eff=0.5, Λ→pπ fraction=0.64이다. 과거 예상 beam 상수는 현재 yield에도 남아 있지만 실제 계산은 hMom2를 쓴다. 출력되는 N_beam 상수값을 입력 ROOT의 적분으로 오해하지 않아야 한다.

## 5. 이 경로로 옮길 때의 구성 제안

```text
sensitivity/
  legacy/polarity_24bin/   # polarity 코드와 해당 param 스냅샷
  legacy/macro_15bin/      # 상위 코드와 해당 param 스냅샷
  data/beam/              # Nbeam CSV/ROOT/convert_csv, 필요하면 과거 beam 자료
  data/mc/                # CSOn, Acceptance_study
  reference/              # total_E72, polar 등 보조 도구와 출처 기록
```

원본을 먼저 보존한 뒤 실행용 코드의 입력·출력 경로와 bin 정의를 통일하는 것을 권장한다. 두 param 폴더는 동명 파일을 합치지 않는다. 핵심 MC 두 개는 약 150 MB이고, polar.cc를 함께 재현할 때만 약 2.70 GB의 CSFlat이 추가로 필요하다. defocus 대용량 자료는 주 파이프라인과 별도다.

현재 폴더에 생성한 `source_inventory.tsv`에는 원본 코드/텍스트 자료의 경로·수정 시각·크기·SHA-256이 있고 `macro_versions.diff`에는 동명 5쌍 및 yield_E72→yield의 실제 차이가 있다. 실제 이관은 아직 수행하지 않았다.
