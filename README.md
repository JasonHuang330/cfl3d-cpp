# cfl3d-cpp

> A C++ translation of NASA's **CFL3D** structured-grid RANS solver.
> NASA CFL3D 结构网格 RANS 求解器的 **C++ 翻译版**。

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![language](https://img.shields.io/badge/C%2B%2B-17-blue.svg)]()

*English* | [中文](#中文说明) | [한국어](#한국어-설명)

---

## English

**cfl3d-cpp** is a C++ port of [CFL3D](https://github.com/nasa/CFL3D) — NASA
Langley's structured-grid, cell-centered, upwind-biased, Reynolds-averaged
Navier-Stokes (RANS) solver (sequential build). The Fortran source was
translated to C++ by an automated Fortran→C++ pipeline, then fixed by hand
until it reproduces the Fortran reference numerically.

It is **~335 translation units** (one `.cpp`/`.h` pair per original routine)
plus a small runtime shim that emulates Fortran array semantics, record-oriented
I/O, and COMMON blocks.

> ⚠️ Independent, unofficial translation. **Not** endorsed by or affiliated with
> NASA. Work in progress — see *Validation* for what is verified.

### Build

```sh
make            # -> ./main   (g++ -O2 -std=c++17)
make clean
```

- Requires a C++17 compiler.
- Requires **libcgns** (for the CGNS grid-I/O path). The Makefile points at
  `/opt/homebrew/{include,lib}`; edit `CGNS_INC` / `CGNS_LIB` if your install
  differs (e.g. `/usr/local`). Cases that don't use CGNS build/run the same.
- macOS only: after a relink you may need `xattr -cr main && codesign -s - -f main`.

### Run

The solver reads the **input deck from stdin** and writes the output files named
in the deck's first "I/O FILES" block:

```sh
./main < your_case.inp
```

### Validation

Verified against the gfortran build of the original CFL3D:

| Case | What it exercises | Result |
|------|-------------------|--------|
| **naca4412** (standard) | steady RANS, single grid | matches Fortran to full precision, `exit 0` |
| **naca4412-xmera** | overset (chimera) | matches Fortran to full precision, `exit 0` |
| **rotorstator** | rotating + overset + dynamic patch + multigrid + unsteady + SA | runs all 1500 steps; residual matches Fortran to 4–6 significant figures |
| **delta** | **CGNS** grid I/O (`icgns=1`), laminar | reads grid from CGNS, runs, residual converges (open item: Cl/Cd print 0) |

Rotorstator was the hardest case; getting it to match uncovered a recurring
translation bug class — *scratch/boundary work-arrays passed to a routine with
the wrong array dimensions* (e.g. `(jdim*kdim, idim*3)` where the callee declares
`(jdim*kdim*idim, 3)`), which corrupts memory via wrong strides. Several such
fixes across the metric / grid-velocity / residual-smoothing paths brought it to
numerical consistency.

### CGNS grid generator

The **delta** case reads its grid from a CGNS file. Generate one from a plot3d
formatted grid plus the deck's BC section:

```sh
g++ -O2 -std=c++17 tools_cgns/gen_delta_cgns.cpp \
    -I/opt/homebrew/include -L/opt/homebrew/lib -lcgns -o tools_cgns/gen_delta_cgns
./tools_cgns/gen_delta_cgns delta.fmt delta_cgns.inp delta.cgns
```

### Tools — grid preprocessors (ronnie, maggie)

`tools/` holds C++ translations of two CFL3D grid-preprocessing utilities. Each
ships its own test case (reference inputs + gfortran goldens) and a
self-contained verify script, so you can reproduce the validation from a fresh
clone with nothing but `g++`:

| Tool | Role | Verified output |
|------|------|-----------------|
| **ronnie** | generalized patched-grid (chimera) connectivity preprocessor | `patch.bin` + `ronnie.out` **byte-exact** vs Fortran |
| **maggie** | overset hole-cutting / interpolation-stencil preprocessor | `ovrlp.bin` **byte-exact** vs Fortran |

```sh
bash tools/ronnie/build_verify.sh
#   compiles, runs the bundled case, byte-compares -> "RESULT: BYTE-EXACT MATCH"

bash tools/maggie/build_verify.sh
#   -> "RESULT: PASS (binary output byte-exact vs Fortran reference)"
```

Each script compiles the port with `g++ -std=c++17`, runs it on the bundled case
in a fresh `run/` sandbox, and byte-compares the result against the gfortran
reference. (maggie's text log also matches, except two "maximum deviation" lines
print the value `0.0` in a different float format — numerically identical.) To
run a tool by hand: `g++ -std=c++17 -O2 tools/ronnie/cpp/*.cpp -o ronnie` then
`./ronnie < tools/ronnie/ref/ronnie.inp` in a directory holding the input
`grid.bin`.

### Layout

Sources are grouped by function under `src/<category>/`, each `foo.cpp` next to
its `foo.h`. The Makefile auto-collects `src/**/*.cpp` and puts every `src/`
subfolder on the include path (sources use flat `#include "foo.h"`), so files can
be added/moved without editing the Makefile. Category is for navigation only and
does not affect the build.

```
src/
  bc/          boundary conditions            io/          restart / plot3d / CGNS / output
  flux/        inviscid + viscous fluxes       force/       forces & diagnostics
  implicit/    approximate-factorization solve aero/        aeroelastic / rigid-body motion
  residual/    residual assembly / update      setup/       setup / sizing / preprocessing
  turbulence/  turbulence models (SA, k-w...)  util/        utilities (parser, complex, alloc)
  metric/      grid metrics & cell volumes     main/        program entry
  grid/        grid read / motion / deform
  overset/     chimera + patched grids
  multigrid/   multigrid transfer
runtime/       FortranArray views, record-aware I/O, COMMON-block structs
tools_cgns/    CGNS grid-file generator (plot3d -> CGNS)
```

### License

Apache License 2.0 — same as upstream CFL3D. See [LICENSE](LICENSE) and
[NOTICE](NOTICE). Original CFL3D © 2001 U.S. Government as represented by the
Administrator of NASA.

---

## 中文说明

**cfl3d-cpp** 是 [CFL3D](https://github.com/nasa/CFL3D) 的 C++ 翻译版 —— CFL3D 是
NASA 兰利研究中心的结构网格、格心格式、迎风偏置的雷诺平均 Navier–Stokes（RANS）
求解器（串行版）。原 Fortran 代码经自动化 Fortran→C++ 流水线翻译，再手工修复，
直到与 Fortran 参考在数值上一致。

项目约 **335 个翻译单元**（每个原例程对应一对 `.cpp`/`.h`），外加一层轻量运行时
垫片，模拟 Fortran 的数组语义、记录式 I/O 和 COMMON 块。

> ⚠️ 独立、非官方翻译，**与 NASA 无关联、未获其背书**。仍在完善中 —— 已验证内容见
> *验证* 一节。

### 构建

```sh
make            # 生成 ./main   (g++ -O2 -std=c++17)
make clean
```

- 需要 C++17 编译器。
- 需要 **libcgns**（用于 CGNS 网格 I/O）。Makefile 默认指向 `/opt/homebrew/{include,lib}`，
  安装路径不同请改 `CGNS_INC` / `CGNS_LIB`（如 `/usr/local`）。不用 CGNS 的算例照常构建运行。
- 仅 macOS：重新链接后可能需要 `xattr -cr main && codesign -s - -f main`。

### 运行

求解器从 **标准输入（stdin）读取输入 deck**，输出文件名取自 deck 第一段 "I/O FILES"：

```sh
./main < 你的算例.inp
```

### 验证

已对照原版 CFL3D 的 gfortran 构建逐一验证：

| 算例 | 考察内容 | 结果 |
|------|----------|------|
| **naca4412**（标准） | 定常 RANS，单网格 | 与 Fortran 全精度一致，`exit 0` |
| **naca4412-xmera** | 重叠网格（chimera） | 与 Fortran 全精度一致，`exit 0` |
| **rotorstator** | 旋转 + 重叠 + 动态拼接 + 多重网格 + 非定常 + SA 湍流 | 跑完全部 1500 步；残差与 Fortran 匹配 4–6 位有效数字 |
| **delta** | **CGNS** 网格 I/O（`icgns=1`），层流 | 从 CGNS 读网格、可运行、残差收敛（遗留项：Cl/Cd 输出为 0） |

rotorstator 是最难的算例；攻克它揭示了一类反复出现的翻译 bug —— *传给例程的临时/
边界工作数组维度算错*（例如声明为 `(jdim*kdim*idim, 3)` 却传成 `(jdim*kdim, idim*3)`），
错误的步长会破坏内存。在度规 / 网格速度 / 残差光顺等路径上修复多处此类 bug 后，
达成数值一致。

### CGNS 网格生成器

**delta** 算例从 CGNS 文件读网格。可用平铺 plot3d 网格 + deck 的 BC 段生成：

```sh
g++ -O2 -std=c++17 tools_cgns/gen_delta_cgns.cpp \
    -I/opt/homebrew/include -L/opt/homebrew/lib -lcgns -o tools_cgns/gen_delta_cgns
./tools_cgns/gen_delta_cgns delta.fmt delta_cgns.inp delta.cgns
```

### 工具 —— 网格预处理器（ronnie、maggie）

`tools/` 收录了 CFL3D 两个网格预处理工具的 C++ 翻译。每个都自带测试算例（参考输入 +
gfortran 金标准输出）和一个自包含的验证脚本，克隆下来只需 `g++` 即可复现验证：

| 工具 | 作用 | 已验证输出 |
|------|------|-----------|
| **ronnie** | 通用拼接网格（chimera）连通性预处理器 | `patch.bin` + `ronnie.out` 与 Fortran **字节完全一致** |
| **maggie** | 重叠网格挖洞 / 插值模板预处理器 | `ovrlp.bin` 与 Fortran **字节完全一致** |

```sh
bash tools/ronnie/build_verify.sh
#   编译、跑算例、字节对比 -> "RESULT: BYTE-EXACT MATCH"

bash tools/maggie/build_verify.sh
#   -> "RESULT: PASS (binary output byte-exact vs Fortran reference)"
```

每个脚本用 `g++ -std=c++17` 编译、在全新的 `run/` 沙箱里跑算例、再与 gfortran 参考做字节对比。
（maggie 的文本日志也一致，只有两行 "maximum deviation" 把数值 `0.0` 用不同浮点格式打印 ——
数值上完全相同。）手动运行:`g++ -std=c++17 -O2 tools/ronnie/cpp/*.cpp -o ronnie`,
然后在放有输入 `grid.bin` 的目录里 `./ronnie < tools/ronnie/ref/ronnie.inp`。

### 目录结构

源码按功能归到 `src/<类别>/`，每个 `foo.cpp` 和它的 `foo.h` 放一起。Makefile 用
`find` 自动收集 `src/**/*.cpp`，并把每个 `src/` 子目录都加进 include path（源码用平铺
的 `#include "foo.h"`），所以加/移文件无需改 Makefile。分类仅为便于导航，不影响编译。

```
src/
  bc/          边界条件            io/          重启 / plot3d / CGNS / 输出
  flux/        无粘 + 粘性通量      force/       气动力与诊断
  implicit/    近似因子分解求解     aero/        气弹 / 刚体运动
  residual/    残差组装 / 更新      setup/       初始化 / 定尺寸 / 预处理
  turbulence/  湍流模型 (SA, k-ω)  util/        工具 (parser, 复数, 内存)
  metric/      网格度规与体积       main/        程序入口
  grid/        网格读取 / 运动 / 变形
  overset/     重叠 + 拼接网格
  multigrid/   多重网格传递
runtime/       FortranArray 视图、记录感知 I/O、COMMON 块结构体
tools_cgns/    CGNS 网格生成器 (plot3d -> CGNS)
```

### 许可证

Apache License 2.0 —— 与上游 CFL3D 相同。见 [LICENSE](LICENSE) 与 [NOTICE](NOTICE)。
原版 CFL3D © 2001 美国政府（由 NASA 局长代表）。

---

## 한국어 설명

**cfl3d-cpp** 는 [CFL3D](https://github.com/nasa/CFL3D) 의 C++ 이식판입니다. CFL3D 는
NASA 랭글리 연구소의 정렬 격자, 셀 중심, 풍상 편향 레이놀즈 평균 나비에–스토크스
(RANS) 솔버(순차 빌드)입니다. 원본 Fortran 코드를 자동 Fortran→C++ 파이프라인으로
번역한 뒤, Fortran 기준 결과와 수치적으로 일치할 때까지 수작업으로 수정했습니다.

약 **335 개의 번역 단위**(원본 루틴마다 `.cpp`/`.h` 한 쌍)와, Fortran 의 배열 의미,
레코드 기반 I/O, COMMON 블록을 흉내 내는 얇은 런타임 계층으로 구성됩니다.

> ⚠️ 독립적·비공식 번역이며 NASA 와 **무관하고 승인받지 않았습니다**. 아직 개선
> 중입니다 — 검증된 범위는 *검증* 절을 참고하세요.

### 빌드

```sh
make            # -> ./main   (g++ -O2 -std=c++17)
make clean
```

- C++17 컴파일러가 필요합니다.
- **libcgns** 가 필요합니다(CGNS 격자 I/O 용). Makefile 은 기본적으로
  `/opt/homebrew/{include,lib}` 를 가리킵니다. 설치 경로가 다르면 `CGNS_INC` /
  `CGNS_LIB` 를 수정하세요(예: `/usr/local`). CGNS 를 쓰지 않는 케이스도 동일하게 빌드·실행됩니다.
- macOS 전용: 재링크 후 `xattr -cr main && codesign -s - -f main` 이 필요할 수 있습니다.

### 실행

솔버는 **입력 덱을 표준 입력(stdin)에서** 읽고, 덱의 첫 "I/O FILES" 블록에 적힌
이름으로 출력 파일을 씁니다:

```sh
./main < your_case.inp
```

### 검증

원본 CFL3D 의 gfortran 빌드와 대조하여 검증했습니다:

| 케이스 | 검증 대상 | 결과 |
|--------|-----------|------|
| **naca4412**(표준) | 정상 RANS, 단일 격자 | Fortran 과 완전 정밀도 일치, `exit 0` |
| **naca4412-xmera** | 중첩 격자(chimera) | Fortran 과 완전 정밀도 일치, `exit 0` |
| **rotorstator** | 회전 + 중첩 + 동적 패치 + 다중 격자 + 비정상 + SA | 1500 스텝 전부 실행; 잔차가 Fortran 과 유효숫자 4–6 자리 일치 |
| **delta** | **CGNS** 격자 I/O(`icgns=1`), 층류 | CGNS 에서 격자 읽기·실행·잔차 수렴(미해결: Cl/Cd 가 0 으로 출력) |

rotorstator 가 가장 어려운 케이스였고, 이를 일치시키는 과정에서 반복되는 번역
버그 유형을 발견했습니다 — *루틴에 넘기는 임시/경계 작업 배열의 차원이 잘못된 경우*
(예: 피호출부가 `(jdim*kdim*idim, 3)` 로 선언했는데 `(jdim*kdim, idim*3)` 로 전달)로,
잘못된 스트라이드가 메모리를 손상시킵니다. 격자 도량·격자 속도·잔차 평활 경로에서
이런 버그 여러 개를 고쳐 수치 일치를 달성했습니다.

### CGNS 격자 생성기

**delta** 케이스는 격자를 CGNS 파일에서 읽습니다. 평문 plot3d 격자와 덱의 BC 섹션으로
생성할 수 있습니다:

```sh
g++ -O2 -std=c++17 tools_cgns/gen_delta_cgns.cpp \
    -I/opt/homebrew/include -L/opt/homebrew/lib -lcgns -o tools_cgns/gen_delta_cgns
./tools_cgns/gen_delta_cgns delta.fmt delta_cgns.inp delta.cgns
```

### 도구 — 격자 전처리기 (ronnie, maggie)

`tools/` 에는 CFL3D 의 격자 전처리 유틸리티 두 개의 C++ 번역이 들어 있습니다. 각각 자체
테스트 케이스(참조 입력 + gfortran 골든 출력)와 자체 완결형 검증 스크립트를 포함하므로,
새로 clone 한 뒤 `g++` 만으로 검증을 재현할 수 있습니다:

| 도구 | 역할 | 검증된 출력 |
|------|------|------------|
| **ronnie** | 일반 패치 격자(chimera) 연결성 전처리기 | `patch.bin` + `ronnie.out` 이 Fortran 과 **바이트 단위 완전 일치** |
| **maggie** | 중첩 격자 홀 커팅 / 보간 스텐실 전처리기 | `ovrlp.bin` 이 Fortran 과 **바이트 단위 완전 일치** |

```sh
bash tools/ronnie/build_verify.sh
#   컴파일·케이스 실행·바이트 비교 -> "RESULT: BYTE-EXACT MATCH"

bash tools/maggie/build_verify.sh
#   -> "RESULT: PASS (binary output byte-exact vs Fortran reference)"
```

각 스크립트는 `g++ -std=c++17` 로 컴파일하고, 새 `run/` 샌드박스에서 케이스를 실행한 뒤
gfortran 참조와 바이트 단위로 비교합니다. (maggie 의 텍스트 로그도 일치하며, 두 개의
"maximum deviation" 줄만 값 `0.0` 을 다른 부동소수점 형식으로 출력합니다 — 수치적으로 동일.)
수동 실행: `g++ -std=c++17 -O2 tools/ronnie/cpp/*.cpp -o ronnie` 후, 입력 `grid.bin` 이
있는 디렉터리에서 `./ronnie < tools/ronnie/ref/ronnie.inp`.

### 디렉터리 구조

소스는 기능별로 `src/<카테고리>/` 아래에 모여 있고, 각 `foo.cpp` 는 짝인 `foo.h` 와
함께 둡니다. Makefile 이 `src/**/*.cpp` 를 자동 수집하고 모든 `src/` 하위 폴더를
include 경로에 추가하므로(소스는 평문 `#include "foo.h"` 사용), Makefile 을 고치지 않고
파일을 추가·이동할 수 있습니다. 카테고리는 탐색용일 뿐 빌드에는 영향을 주지 않습니다.

```
src/
  bc/          경계 조건            io/          재시작 / plot3d / CGNS / 출력
  flux/        비점성 + 점성 플럭스  force/       공력 하중 및 진단
  implicit/    근사 인수분해 솔버    aero/        공탄성 / 강체 운동
  residual/    잔차 조립 / 갱신      setup/       설정 / 크기 산정 / 전처리
  turbulence/  난류 모델 (SA, k-ω)  util/        유틸리티 (parser, 복소수, 할당)
  metric/      격자 도량과 셀 부피   main/        프로그램 진입점
  grid/        격자 읽기 / 이동 / 변형
  overset/     중첩 + 패치 격자
  multigrid/   다중 격자 전달
runtime/       FortranArray 뷰, 레코드 인식 I/O, COMMON 블록 구조체
tools_cgns/    CGNS 격자 생성기 (plot3d -> CGNS)
```

### 라이선스

Apache License 2.0 —— 상류 CFL3D 와 동일합니다. [LICENSE](LICENSE) 및 [NOTICE](NOTICE)
참고. 원본 CFL3D © 2001 미국 정부(NASA 국장 대표).
