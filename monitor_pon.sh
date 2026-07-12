#!/bin/bash
#
# Copyright (c) 2026, Po-Yu Lin (林伯昱)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# THIRD-PARTY TOOLS NOTICE:
# Users of main_research.nf must comply with:
#   - Manta: PolyForm Strict License 1.0.0 (non-commercial only)
#   - ExpansionHunter: PolyForm Strict License 1.0.0 (non-commercial only)
# See README.md and LICENSE for details.

# =====================================================================
# DGX-2 PON Monitor | 用法：bash monitor_pon.sh
#
# 可用環境變數覆寫（方便沿用到 case pipeline 或其他機器，避免再度寫死路徑）：
#   PON_WORK_DIR      啟動目錄（.nextflow.log 所在），預設 /raid/DGM/pon_work
#   GPU_LOCK_DIR      GPU lock 目錄，需與 gpu_lock.sh 一致，預設 /raid/DGM/gpu_locks
#   MONITOR_GPU_IDS   監看的 GPU id（空白分隔），預設 "10 11 12 13 14 15"
#   MONITOR_INTERVAL  更新間隔秒數，預設 30
#   MONITOR_NXF_MATCH pgrep 比對的 nextflow 進程樣式，預設 "nextflow.*main_pon"
# =====================================================================

# ── 可覆寫設定 ────────────────────────────────────────────────
WORK_DIR="${PON_WORK_DIR:-/raid/DGM/pon_work}"

# ★ 修正：GPU lock 目錄必須與 gpu_lock.sh / gpu_unlock.sh 一致（/raid/DGM/gpu_locks）。
#    先前寫死成 /tmp/nxf_gpu_locks，永遠找不到 lock，🔒 一直顯示未鎖定（假象：卡明明忙）。
GPU_LOCK_DIR="${GPU_LOCK_DIR:-/raid/DGM/gpu_locks}"

INTERVAL="${MONITOR_INTERVAL:-30}"
NXF_MATCH="${MONITOR_NXF_MATCH:-nextflow.*main_pon}"
LOG="${WORK_DIR}/.nextflow.log"

if [ -n "${MONITOR_GPU_IDS:-}" ]; then
    read -ra GPU_IDS <<< "${MONITOR_GPU_IDS}"
else
    GPU_IDS=(10 11 12 13 14 15)
fi

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; DIM='\033[2m'; NC='\033[0m'

PROCS=("FASTP" "PARABRICKS_FQ2BAM" "PREP_GATK_INTERVALS" "PREP_CNVKIT_BEDS"
       "COLLECT_GATK_COUNTS" "COLLECT_CNVKIT_COV" "CNVKIT_REFERENCE"
       "FILTER_INTERVALS" "PLOIDY_COHORT" "SCATTER_INTERVALS" "GCNV_COHORT")

# 短名對照（控制欄寬）
declare -A SHORT=(
    [FASTP]="FASTP"
    [PARABRICKS_FQ2BAM]="FQ2BAM"
    [PREP_GATK_INTERVALS]="PREP_GATK"
    [PREP_CNVKIT_BEDS]="PREP_CNVKIT"
    [COLLECT_GATK_COUNTS]="COLL_GATK"
    [COLLECT_CNVKIT_COV]="COLL_CNVKIT"
    [CNVKIT_REFERENCE]="CNVKIT_REF"
    [FILTER_INTERVALS]="FILTER_INT"
    [PLOIDY_COHORT]="PLOIDY"
    [SCATTER_INTERVALS]="SCATTER"
    [GCNV_COHORT]="GCNV"
)

# 從 .nextflow.log 解析單一 process 的計數，回傳：submitted completed running failed cached
#   run    = submitted - completed - failed（實際執行中）
#   cached = -resume 時被快取略過的任務（先前已完成，不會再出現 status 行）
# 註：Nextflow log 沒有可靠的「排隊等待」訊號（原本的 Wait 欄公式恆為 0），故改用 Cache。
get_stats() {
    local p=$1
    [ ! -f "$LOG" ] && echo "0 0 0 0 0" && return
    local sub cac com fai
    sub=$(grep -c "Submitted process > ${p}"          "$LOG" 2>/dev/null)
    cac=$(grep -c "Cached process > ${p}"             "$LOG" 2>/dev/null)
    # 用 [; ] 當分隔，同時相容「有樣本名」與「無樣本名」兩種 log 格式
    com=$(grep -c "name: ${p}[; ].*status: COMPLETED" "$LOG" 2>/dev/null)
    fai=$(grep -c "name: ${p}[; ].*status: FAILED"    "$LOG" 2>/dev/null)
    sub=${sub//[^0-9]/}; sub=${sub:-0}
    cac=${cac//[^0-9]/}; cac=${cac:-0}
    com=${com//[^0-9]/}; com=${com:-0}
    fai=${fai//[^0-9]/}; fai=${fai:-0}
    local run=$(( sub - com - fai ))
    [ $run -lt 0 ] && run=0
    echo "$sub $com $run $fai $cac"
}

bar() {
    local cur=${1:-0} tot=${2:-0} w=${3:-10}
    cur=${cur//[^0-9]/}; cur=${cur:-0}
    tot=${tot//[^0-9]/}; tot=${tot:-0}
    local f=0; [ $tot -gt 0 ] && f=$(( cur * w / tot ))
    [ $f -gt $w ] && f=$w
    local e=$(( w - f )) b=""
    for((i=0;i<f;i++)); do b+="█"; done
    for((i=0;i<e;i++)); do b+="░"; done
    echo "$b"
}

CYCLE=0
WS="?"

while true; do
    clear

    # ── Header ──────────────────────────────────────────────────
    NXF_PID=$(pgrep -f "${NXF_MATCH}" 2>/dev/null | head -1)
    [ -n "$NXF_PID" ] && NXF_ST="${GREEN}●RUN${NC}(${NXF_PID})" || NXF_ST="${RED}●STOP${NC}"
    printf "${BOLD}[ DGX-2 PON Monitor ]${NC} $(date '+%m-%d %H:%M:%S')  Nextflow: "
    echo -e "${NXF_ST}"

    # ── Process 進度（壓縮成 1 行/process）──────────────────────
    echo -e "${DIM}$(printf '─%.0s' {1..90})${NC}"
    # 欄位：Process / 完成(含快取) / 執行中 / 快取略過 / 失敗 / 進度
    printf "  ${BOLD}%-12s %5s %5s %5s %5s  %-10s${NC}\n" \
        "Process" "Done" "Run" "Cache" "Fail" "Progress"
    echo -e "  ${DIM}$(printf '─%.0s' {1..86})${NC}"

    for p in "${PROCS[@]}"; do
        read -r sub com run fai cac <<< "$(get_stats "$p")"
        sub=${sub:-0}; com=${com:-0}; run=${run:-0}; fai=${fai:-0}; cac=${cac:-0}
        # Done = 已完成（新跑完 + 快取略過）；Total = 已提交 + 快取
        total=$(( sub + cac ))
        fin=$(( com + cac ))

        if   [ "$fai" -gt 0 ]; then sc=$RED
        elif [ "$run" -gt 0 ]; then sc=$CYAN
        elif [ "$total" -gt 0 ] && [ "$fin" -eq "$total" ]; then sc=$GREEN
        else sc=$DIM; fi

        b=$(bar "$fin" "$total" 10)
        [ $total -gt 0 ] && prog="${fin}/${total}" || prog="-"

        printf "  ${sc}%-12s${NC} %5s ${CYAN}%5s${NC} ${DIM}%5s${NC} " \
            "${SHORT[$p]}" "$fin" "$run" "$cac"
        [ "$fai" -gt 0 ] && printf "${RED}%5s${NC}  " "$fai" || printf "%5s  " "$fai"

        if [ "$total" -gt 0 ] && [ "$fin" -eq "$total" ]; then
            echo -e "${GREEN}${b}${NC} ${prog}"
        elif [ "$run" -gt 0 ]; then
            echo -e "${CYAN}${b}${NC} ${prog}"
        else
            echo -e "${DIM}${b}${NC} ${prog}"
        fi
    done

    # ── 系統資源（單行）─────────────────────────────────────────
    echo -e "${DIM}$(printf '─%.0s' {1..90})${NC}"

    # CPU
    CPU_USAGE=0
    if command -v top >/dev/null 2>&1; then
        CPU_USAGE=$(top -bn1 | awk '/^%Cpu/{print 100-$8}' | cut -d. -f1 2>/dev/null)
    fi
    CPU_USAGE=${CPU_USAGE//[^0-9]/}; CPU_USAGE=${CPU_USAGE:-0}
    # MEM
    read -r MT MU <<< "$(free -g | awk '/^Mem:/{print $2,$3}')"
    MT=${MT:-1}; MU=${MU:-0}; [ "$MT" -eq 0 ] && MT=1; MP=$(( MU * 100 / MT ))
    printf "  CPU $(bar $CPU_USAGE 100 10) ${CPU_USAGE}%%  |  MEM $(bar $MU $MT 10) ${MU}/${MT}GB (${MP}%%)\n"

    # ── GPU（3 個一排）──────────────────────────────────────────
    # ★ 優化：一次 nvidia-smi 查全部 GPU（原本每張卡各叫一次，每次更新 6 次 fork）
    echo -e "${DIM}$(printf '─%.0s' {1..90})${NC}"
    unset GVU GVT GUT; declare -A GVU GVT GUT
    if command -v nvidia-smi >/dev/null 2>&1; then
        while IFS=',' read -r gidx gused gtot gutil; do
            gidx=${gidx//[^0-9]/}
            [ -z "$gidx" ] && continue
            GVU[$gidx]=${gused//[^0-9]/}; GVT[$gidx]=${gtot//[^0-9]/}; GUT[$gidx]=${gutil//[^0-9]/}
        done < <(nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu \
                    --format=csv,noheader,nounits -i "$(IFS=,; echo "${GPU_IDS[*]}")" 2>/dev/null)
    fi

    COUNT=0
    for GPU_ID in "${GPU_IDS[@]}"; do
        VU=${GVU[$GPU_ID]:-0}; VT=${GVT[$GPU_ID]:-1}; UT=${GUT[$GPU_ID]:-0}
        VU=${VU:-0}; VT=${VT:-1}; UT=${UT:-0}; [ "$VT" -eq 0 ] && VT=1
        VUG=$(( VU / 1024 )); VTG=$(( VT / 1024 ))
        # 🔒 = 這張卡目前被 gpu_lock.sh 鎖定（有 gpu_<id>.lock）
        if [ -f "${GPU_LOCK_DIR}/gpu_${GPU_ID}.lock" ]; then
            LK="${GREEN}🔒${NC}"; vc=$CYAN
        else
            LK="${DIM}--${NC}"; vc=$DIM
        fi
        printf "  GPU%-2s ${vc}%3s/%-3sGB${NC} %3s%% ${LK}" "$GPU_ID" "$VUG" "$VTG" "$UT"
        COUNT=$(( COUNT + 1 ))
        [ $(( COUNT % 3 )) -eq 0 ] && echo "" || printf "   "
    done
    # 若最後一排未滿 3 張，補一個換行
    [ $(( COUNT % 3 )) -ne 0 ] && echo ""

    # ── 磁碟（單行）─────────────────────────────────────────────
    # ★ 優化：du 會 stat 整棵 work 樹，改成每 10 圈（約 5 分鐘）才重算一次，
    #    避免每 30s 對 /raid SSD 造成 stat 風暴、與 pipeline 的 I/O 互搶。
    #    df 是即時的、幾乎零成本，所以每圈都更新。
    echo -e "${DIM}$(printf '─%.0s' {1..90})${NC}"
    read -r DU DT DA DP <<< "$(df -h /raid | awk 'NR==2{gsub(/%/,"",$5); print $3,$2,$4,$5}')"
    DP=${DP:-0}
    if [ $(( CYCLE % 10 )) -eq 0 ]; then
        WS=$(du -sh "${WORK_DIR}" 2>/dev/null | cut -f1)
        WS=${WS:-?}
    fi
    echo -e "  /raid $(bar $DP 100 10) ${DU}/${DT} 剩${DA} (${DP}%)  work:${CYAN}${WS}${NC}  ${DIM}[Ctrl+C 離開]${NC}"

    CYCLE=$(( CYCLE + 1 ))
    sleep "${INTERVAL}"
done
