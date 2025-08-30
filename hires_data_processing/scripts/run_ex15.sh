#!/usr/bin/env bash
# Unified runner for schicluster pipeline (per-stage logging)
# 默认: 处理 E70，分辨率 100000，chr1=1 pos1=2 chr2=3 pos2=4
# 仅支持 micromamba 环境；若不存在则报错退出

# -------------------------------
# Defaults (可根据需要修改)
# -------------------------------
STAGES="EX15"                       # 默认处理阶段 E70
BASE_DIR=""                        # 可用 -b 指定；为空则使用脚本相对目录的上级
ENV_NAME="schicluster"             # 可用 -n 指定
PYTHON_BIN="python3"               # 可用 -p 指定
SCRIPT_REL="scripts/run_schicluster_pipeline.py"  # Python 脚本相对 BASE_DIR 的路径
RESOLUTION="100000"                # 默认分辨率 100000，可用 -r 覆盖
BATCH_SIZE=""                      # 可用 -B 指定
CPU_PER_JOB=""                     # 可用 -c 指定
CHR1="1"                           # 默认 chr1 列索引
POS1="2"                           # 默认 pos1 列索引
CHR2="3"                           # 默认 chr2 列索引
POS2="4"                           # 默认 pos2 列索引
DRY_RUN=0                          # 可用 -d 开启
EXTRA_ARGS=()                      # 透传给 Python 的其它参数（放在 -- 之后）

# -------------------------------
# Parse args
# -------------------------------
print_help() {
  cat <<'EOF'
Options:
  -s "E70 E75 EX05"   空格分隔的阶段列表（默认：E70）
  -b /abs/base_dir    项目根目录；不填则用脚本所在目录的上级
  -n env_name         micromamba 环境名（默认：schicluster）
  -p /path/python     Python 解释器（默认：python3）
  -r 1000000          分辨率；传给 Python 的 --resolution（默认：100000）
  -B 2000             批量大小；传给 Python 的 --batch_size
  -c 30               CPU 数；传给 Python 的 --cpu_per_job
  -d                  dry-run（只打印命令不执行）
  --                  之后的参数全部原样透传给 Python 脚本
EOF
}

PASSTHROUGH=0
while (( "$#" )); do
  if [[ $PASSTHROUGH -eq 1 ]]; then
    EXTRA_ARGS+=("$1"); shift; continue
  fi
  case "${1:-}" in
    -s) STAGES="$2"; shift 2 ;;
    -b) BASE_DIR="$2"; shift 2 ;;
    -n) ENV_NAME="$2"; shift 2 ;;
    -p) PYTHON_BIN="$2"; shift 2 ;;
    -r) RESOLUTION="$2"; shift 2 ;;
    -B) BATCH_SIZE="$2"; shift 2 ;;
    -c) CPU_PER_JOB="$2"; shift 2 ;;
    -d) DRY_RUN=1; shift ;;
    -h|--help) print_help; exit 0 ;;
    --) PASSTHROUGH=1; shift ;;
    *) echo "Unknown option: $1"; print_help; exit 2 ;;
  esac
done

# -------------------------------
# Resolve BASE_DIR & paths
# -------------------------------
if [[ -z "$BASE_DIR" ]]; then
  SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
  BASE_DIR="$( dirname "$SCRIPT_DIR" )"
fi

PROJECT_DIR="$BASE_DIR"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

PY_SCRIPT="$PROJECT_DIR/$SCRIPT_REL"
if [[ ! -f "$PY_SCRIPT" ]]; then
  ALT="$PROJECT_DIR/run_schicluster_pipeline.py"
  if [[ -f "$ALT" ]]; then
    PY_SCRIPT="$ALT"
  else
    echo "ERROR: Python script not found: $PY_SCRIPT (or $ALT)"
    exit 1
  fi
fi

# -------------------------------
# Activate environment (micromamba ONLY)
# -------------------------------
# Activate environment and run processing
eval "$(micromamba shell hook --shell=bash)"
micromamba activate schicluster

# -------------------------------
# Build common args for Python
# -------------------------------
COMMON_ARGS=( )
COMMON_ARGS+=( "--base_dir" "$PROJECT_DIR" )
[[ -n "$RESOLUTION"  ]] && COMMON_ARGS+=( "--resolution" "$RESOLUTION" )
[[ -n "$BATCH_SIZE"  ]] && COMMON_ARGS+=( "--batch_size" "$BATCH_SIZE" )
[[ -n "$CPU_PER_JOB" ]] && COMMON_ARGS+=( "--cpu_per_job" "$CPU_PER_JOB" )
COMMON_ARGS+=( "--chr1" "$CHR1" "--pos1" "$POS1" "--chr2" "$CHR2" "--pos2" "$POS2" )
(( DRY_RUN == 1 )) && COMMON_ARGS+=( "--dry_run" )

# -------------------------------
# Run per stage
# -------------------------------
TS="$(date +%Y%m%d_%H%M%S)"
echo "Starting processing at $(date)"
echo "Stages: $STAGES"
echo "Logs:   $LOG_DIR"
echo "Base:   $PROJECT_DIR"
echo "Script: $PY_SCRIPT"
echo "Default resolution: $RESOLUTION"
echo "Default columns: chr1=$CHR1 pos1=$POS1 chr2=$CHR2 pos2=$POS2"
echo

EXIT_SUM=0
for STAGE in $STAGES; do
  LOG_FILE="$LOG_DIR/processing_${STAGE}_${TS}.log"
  echo "----- Stage $STAGE START $(date) -----" | tee -a "$LOG_FILE"

  set +e
  CMD=( "$PYTHON_BIN" "$PY_SCRIPT" "${COMMON_ARGS[@]}" "--specific_stage" "$STAGE" "${EXTRA_ARGS[@]}" )
  echo "CMD: ${CMD[*]}" | tee -a "$LOG_FILE"

  "${CMD[@]}" 2>&1 | tee -a "$LOG_FILE"
  STAGE_EXIT=${PIPESTATUS[0]}

  if [[ $STAGE_EXIT -eq 0 ]]; then
    echo "Stage $STAGE SUCCESS $(date)" | tee -a "$LOG_FILE"
    touch "$LOG_DIR/completed_${STAGE}.flag"
  else
    echo "Stage $STAGE FAILED (code=$STAGE_EXIT) $(date)" | tee -a "$LOG_FILE"
    touch "$LOG_DIR/failed_${STAGE}.flag"
    EXIT_SUM=$STAGE_EXIT
  fi
  echo "----- Stage $STAGE END -----" | tee -a "$LOG_FILE"
  echo
  set -e
done

# -------------------------------
# Deactivate & exit
# -------------------------------
micromamba deactivate
exit $EXIT_CODE