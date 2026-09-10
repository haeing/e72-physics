# E72 sensitivity — 실제 실험 빔 입력

`Nbeam_251123.root`가 실제 실험에서 사용한 빔 수 자료라는 사용자 확인에 따라, 해당 파일의 `hMom2`를 입력으로 최신 상위 macro의 **15 운동량 bin** 계산을 이관했다. 원본은 그대로 두고 코드와 필수 데이터를 이 디렉터리에 복사했다.

## 실행

ROOT가 PATH에 있는 환경에서:

```bash
bash run_all.sh
```

다른 디렉터리에서 스크립트 절대경로로 실행해도 된다. 실행 시 `results/`의 ROOT 결과와 합본 PDF를 갱신한다. TXT 파일은 읽거나 생성하지 않는다. 각 단계는 전역변수 충돌을 피하도록 별도 ROOT 프로세스로 실행한다.

개별 그림을 대화형으로 보려면 이 디렉터리에서 `root -l macro/polar_th.cc`와 같이 실행한다. 하위 단계에는 앞 단계의 ROOT 결과가 필요하다.

## 유지한 계산과 결과

| 실행 순서 | 코드 | 입력/결과 |
|---|---|---|
| 1 | `macro/yield.cc` | 실제 빔 `data/beam/Nbeam_251123.root:hMom2` → `results/yield.root`의 `entry_mom`, `luminosity_mom` 그래프 |
| 2 | `macro/acceptance.cc` | `data/mc/7201_7202_Acceptance_study.root` → `results/acceptance.root`의 운동량별 acceptance 그래프 |
| 3 | `macro/diff_cross.cc` | `data/mc/7201_7202_CSOn.root` + entry_mom → `results/diff_cross.root`의 운동량별 yield 그래프 (각도별 예상 yield) |
| 4 | `macro/polar_th.cc` | 각도별 yield → Λ 편극 모델 비교 및 통계 sensitivity |
| 5 | `macro/crystalball.cc` | yield, luminosity, acceptance + 문헌 배열 → differential cross section 및 예상 통계오차 |

전체 그림은 실행 순서(yield → acceptance → diff_cross → polar_th → crystalball)대로 **`results/sensitivity.pdf` 하나의 multipage PDF**에 저장된다. 각 canvas가 한 페이지이며, PDF 북마크에 매크로와 canvas 이름을 표시한다. 매 실행 시 이 PDF를 갱신한다. 편집 가능한 canvas는 `results/<매크로>.root`, 로그는 `results/<매크로>.log`에 저장된다. 기존 ROOT 결과만 다시 PDF로 묶으려면 `root -l -b -q macro/export_pdf.C`를 실행한다.

`macro/total_E72.cc`도 문헌 total cross section 참고용으로 복사했다. 전체 실행에는 포함하지 않았다.

## ROOT 그래프 사용

canvas 외에 그래프를 **독립된 `TGraphErrors` 객체**로 저장한다. 운동량 디렉터리의 `p724`는 중심 724 MeV/c, 반폭 1 MeV/c이다. 15개 운동량 모두 동일한 구조다.

| ROOT 파일 | 객체 경로 | 내용 |
|---|---|---|
| `results/yield.root` | `entry_mom`, `luminosity_mom` | 운동량별 yield / luminosity |
| `results/yield.root` | `yield_spectrum`, `luminosity_spectrum`, `beam_momentum` | 전체 스펙트럼 그래프와 입력 빔 히스토그램 |
| `results/acceptance.root` | `momentum/p724/acceptance` | cosθ별 acceptance와 기존 오차 |
| `results/diff_cross.root` | `momentum/p724/yield` | cosθ별 예상 사건 수 |
| `results/crystalball.root` | `momentum/p724/differential_cs` | cosθ별 단면적과 예상 통계오차 |
| `results/crystalball.root` | `momentum/p724/crystal_ball`, `angular_fit` | 문헌 그래프 및 fit 함수 |
| `results/polar_th.root` | `momentum/p724/polarization_uncertainty` | cosθ별 편극 통계오차 (y값이 δP) |
| `results/polar_th.root` | `polarization_uncertainty_momentum` | 운동량별 각도 합산 편극 통계오차 |
| `results/polar_th.root` | `model_comparison/model735_sum724_750/` 및 `model765_bin766/` | Crystal Ball, D03/P03 모델 그래프와 spline |

예를 들어 ROOT에서:

```cpp
auto f = TFile::Open("results/acceptance.root");
auto g = f->Get<TGraphErrors>("momentum/p724/acceptance");
g->Draw("AP");
```

acceptance/yield/단면적/편극오차 각도 그래프는 같은 디렉터리의 `h_<이름>`에 TH1D로도 저장한다. acceptance의 `mc_pass`, `mc_total`은 효율 계산에 사용한 원래 MC 사건 수 히스토그램이다. `diff_cross`에도 `mc_pass`를 저장한다. 각 단계의 canvas는 기존처럼 ROOT 최상위에 남는다.

단계 간 입력은 ROOT 그래프를 직접 읽으며 운동량·각도 좌표 및 점 수를 검사한다. 기존 TXT 출력은 `legacy/last_txt_results/`에 보관했다. 과거 각도별 yield 그래프의 y오차 0 설정은 유지했으며, 편극과 단면적 통계오차는 N에서 별도로 계산한다. 운동량별 yield 그래프에는 sqrt(N) 오차를 저장한다.

## 입력 및 재현성

- 실제 빔 ROOT/CSV 및 변환 코드는 `data/beam/`에 있다. ROOT `hMom2`의 bin은 575~1003 MeV/c, 폭 2 MeV/c이며 적분은 약 1.94083×10⁹이다.
- 필수 MC 두 파일은 `data/mc/`에 실제 복사했다. 외부 경로를 가리키는 링크가 아니다.
- 현재 계산은 724,726,728,730,732,734,738,742,746,750,754,758,762,766,770 MeV/c 중심, 각각 ±1 MeV/c, 각도 15-bin이다.
- 원본 상위 매크로 및 param은 `legacy/macro_15bin/`, polarity 계열은 `legacy/polarity_24bin/`에 수정 없이 보존했다. legacy는 기록용이며 외부 상대경로와 과거 코드 문제가 남아 있다.
- 예상 가우시안 빔의 `yield_E72.cc`는 legacy에만 두었다. 기본 실행은 실제 빔 입력만 사용한다.
- 직접 입자별 편극을 계산하던 `polar.cc`도 legacy에 보존했다. 별도 2.70 GB CSFlat MC는 복사하지 않았으며 현재 편극 sensitivity 계산에는 필요하지 않다.

## 이관 시 변경한 부분

입력 경로를 새 데이터 위치에 맞추고, 15-bin의 운동량 빈틈 및 최상단 경계에서 잘못된 배열 인덱스를 사용하지 않도록 검사했다. yield가 입력 빔 히스토그램을 잘못된 bin 번호로 덮어쓰던 코드를 제거하고 출력 N_beam을 실제 히스토그램 적분으로 변경했다.

crystalball의 보조 pseudo-global fit 실행은 비활성화했다. 기존 함수는 x[1]을 읽지만 호출하는 1차원 TGraph fit은 해당 좌표를 제공하지 않는다. 개별 각도분포 fit 및 differential cross section sensitivity 계산은 유지했다. 원래 코드는 legacy에서 확인할 수 있다.

물리 설정 및 기존 계산식은 유지했다: E_eff=0.8, trigger_eff=0.5, Λ→pπ fraction=0.64, alpha=0.65. 기존 yield의 bin 폭 정규화, crystalball의 Legendre2 정수 나눗셈, acceptance 표시 오차식은 이번 이관에서 변경하지 않았다. 따라서 결과는 **최신 기존 분석을 이관한 결과**이며 물리 계산식 검증까지 끝났다는 의미는 아니다. 확인 사항은 `MACRO_HISTORY.md`의 4절에 있다.

특히 편극의 735 MeV/c 모델 비교 오차에는 첫 10개 운동량 bin(중심 724~750)의 합산 통계를 사용한다. 765 MeV/c 모델에는 중심 766 bin 통계를 사용한다. 이는 최신 상위 코드의 동작을 그대로 유지한 것이다.

## 검증

ROOT 6.32.08에서 다섯 단계의 batch 실행 완료. 15페이지 multipage PDF 하나와 5개 ROOT 결과를 생성한다. 운동량별 TGraphErrors 60개, 유한한 값, acceptance 범위, 각도별 yield 합, 편극 및 단면적 오차식과 저장값의 일치를 검사했다. 이전 TXT 결과와의 비교도 통과했다(상대 허용오차 1e-5; ROOT는 TXT 소수점 반올림 없이 double 값을 전달). `root -l -b -q macro/validate_root.C`로 다시 검사할 수 있다. 복사한 네 입력 데이터는 원본과 동일함을 확인했다. 상세 기록과 SHA-256은 `results/validation.json`에 있다.

과거 계보는 `MACRO_HISTORY.md`, 원본 코드 비교는 `macro_versions.diff`, 조사 당시 원본 목록은 `source_inventory.tsv`를 참고한다. 계보 문서의 '이관 미수행' 표현은 최초 조사 시점의 기록이며, 현재 이관 상태는 이 README가 기준이다.
