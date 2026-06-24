#!/bin/bash

set -u

RESULTS_FILE="PM2DGEMMResultsCANNONS.csv"
RAW_LOG_DIR="dgemmRawLogsCANNONS"

DEFAULT_EXECUTABLE="./pm2dgemmCANNONS"

CHECKSUM_ABS_TOLERANCE="0.000001"
CHECKSUM_REL_TOLERANCE="0.000000001"

declare -A CHECKSUM_BASELINES

HOSTS=(
    hermes1.uwb.edu
    hermes2.uwb.edu
    hermes3.uwb.edu
    hermes4.uwb.edu
    hermes5.uwb.edu
    hermes6.uwb.edu
    hermes7.uwb.edu
    hermes8.uwb.edu
    hermes9.uwb.edu
    hermes10.uwb.edu
    hermes11.uwb.edu
    hermes12.uwb.edu
    hermes13.uwb.edu
    hermes14.uwb.edu
    hermes15.uwb.edu
    hermes16.uwb.edu
    hermes17.uwb.edu
    hermes18.uwb.edu
    hermes19.uwb.edu
    hermes20.uwb.edu
    hermes21.uwb.edu
    hermes22.uwb.edu
    hermes23.uwb.edu
    hermes24.uwb.edu
)

NODE_COUNTS=(1 2 4 8 16 24)
PPN_VALUES=(1 2 4 6)

RUNS=4
SEED=1

SCENARIO_NAMES=(S1Medium S2Large)
SCENARIO_M=(8192 11520)
SCENARIO_N=(8192 11520)
SCENARIO_K=(8192 11520)

read -r -p "CANNONS executable [${DEFAULT_EXECUTABLE}]: " EXECUTABLE
EXECUTABLE="${EXECUTABLE:-$DEFAULT_EXECUTABLE}"

mkdir -p "$RAW_LOG_DIR"
: > "$RESULTS_FILE"

printf 'algorithm,scenario,M,N,K,nodes,ppn,ranks,Px,Py,runsCompleted,avgComputeSeconds,avgWaitSeconds,avgFillSeconds,avgTotalSeconds,checksum,status,logfile\n' \
    >> "$RESULTS_FILE"

buildNodeList() {
    local nodeCount="$1"
    local ppn="$2"
    local nodeList=""
    local i
    local j

    if (( nodeCount > ${#HOSTS[@]} )); then
        echo "ERROR: requested $nodeCount nodes, but only ${#HOSTS[@]} hosts are defined" >&2
        return 1
    fi

    for (( i = 0; i < nodeCount; i++ )); do
        for (( j = 0; j < ppn; j++ )); do
            if [[ -n "$nodeList" ]]; then
                nodeList+=","
            fi
            nodeList+="${HOSTS[$i]}"
        done
    done

    printf '%s\n' "$nodeList"
}

pickBestDecomp() {
    local ranks="$1"
    local bestPx="$ranks"
    local bestPy=1
    local bestDiff="$ranks"
    local py
    local px
    local diff

    for (( py = 1; py * py <= ranks; py++ )); do
        if (( ranks % py != 0 )); then
            continue
        fi

        px=$(( ranks / py ))
        diff=$(( px - py ))

        if (( diff < bestDiff )); then
            bestDiff="$diff"
            bestPx="$px"
            bestPy="$py"
        fi
    done

    printf '%d %d\n' "$bestPx" "$bestPy"
}

pickCannonsDecomp() {
    local ranks="$1"
    local p=1

    while (( p * p < ranks )); do
        ((p++))
    done

    if (( p * p == ranks )); then
        printf '%d %d\n' "$p" "$p"
        return 0
    fi

    return 1
}

safeLogName() {
    local algorithm="$1"
    local scenarioName="$2"
    local matrixSize="$3"
    local nodeCount="$4"
    local ppn="$5"
    local totalRanks="$6"
    local px="$7"
    local py="$8"

    printf '%s/%s_%s_N%s_nodes%s_ppn%s_ranks%s_Px%s_Py%s.log\n' \
        "$RAW_LOG_DIR" \
        "$algorithm" \
        "$scenarioName" \
        "$matrixSize" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks" \
        "$px" \
        "$py"
}

parseMetrics() {
    local output="$1"

    printf '%s\n' "$output" | awk '
        /^Compute Time:/ {
            value = $3
            sub(/s$/, "", value)
            compute += value
            computeCount++
        }
        /^Broadcast Wait Time:/ {
            value = $4
            sub(/s$/, "", value)
            wait += value
            waitCount++
        }
        /^Shift Wait Time:/ {
            value = $4
            sub(/s$/, "", value)
            wait += value
            waitCount++
        }
        /^Fill Time:/ {
            value = $3
            sub(/s$/, "", value)
            fill += value
            fillCount++
        }
        /^Total MPI Timed Phase:/ {
            value = $5
            sub(/s$/, "", value)
            total += value
            totalCount++
        }
        END {
            if (computeCount > 0 && waitCount > 0 && fillCount > 0 && totalCount > 0) {
                printf "%.6f %.6f %.6f %.6f %d\n", \
                    compute / computeCount, \
                    wait / waitCount, \
                    fill / fillCount, \
                    total / totalCount, \
                    computeCount
            }
        }
    '
}

parseChecksum() {
    local output="$1"

    printf '%s\n' "$output" | awk '
        tolower($0) ~ /checksum/ {
            for (i = 1; i <= NF; i++) {
                value = $i
                gsub(/^[^0-9+.-]*/, "", value)
                gsub(/[^0-9eE+.-]*$/, "", value)

                if (value ~ /^[-+]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][-+]?[0-9]+)?$/) {
                    print value
                    exit
                }
            }
        }
    '
}

compareChecksums() {
    local actual="$1"
    local expected="$2"

    awk \
        -v actual="$actual" \
        -v expected="$expected" \
        -v absTol="$CHECKSUM_ABS_TOLERANCE" \
        -v relTol="$CHECKSUM_REL_TOLERANCE" '
        BEGIN {
            diff = actual - expected
            if (diff < 0) {
                diff = -diff
            }

            threshold = absTol + relTol * (expected < 0 ? -expected : expected)

            if (diff <= threshold) {
                print "MATCH"
            } else {
                print "NO_MATCH"
            }
        }
    '
}

setChecksumCsvValue() {
    local algorithm="$1"
    local scenarioName="$2"
    local nodeCount="$3"
    local ppn="$4"
    local checksum="$5"
    local baselineKey="${algorithm}_${scenarioName}"
    local baseline

    if [[ -z "$checksum" ]]; then
        CHECKSUM_VALUE="NO_CHECKSUM"
        return 0
    fi

    if (( nodeCount == 1 && ppn == 1 )); then
        CHECKSUM_BASELINES["$baselineKey"]="$checksum"
        CHECKSUM_VALUE="$checksum"
        return 0
    fi

    baseline="${CHECKSUM_BASELINES[$baselineKey]:-}"

    if [[ -z "$baseline" ]]; then
        CHECKSUM_VALUE="NO_BASELINE"
        return 0
    fi

    CHECKSUM_VALUE="$(compareChecksums "$checksum" "$baseline")"
}

writeResultRow() {
    local algorithm="$1"
    local scenarioName="$2"
    local m="$3"
    local n="$4"
    local k="$5"
    local nodeCount="$6"
    local ppn="$7"
    local totalRanks="$8"
    local px="$9"
    local py="${10}"
    local runsCompleted="${11}"
    local avgComputeSeconds="${12}"
    local avgWaitSeconds="${13}"
    local avgFillSeconds="${14}"
    local avgTotalSeconds="${15}"
    local checksumValue="${16}"
    local status="${17}"
    local logfile="${18}"

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$algorithm" \
        "$scenarioName" \
        "$m" \
        "$n" \
        "$k" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks" \
        "$px" \
        "$py" \
        "$runsCompleted" \
        "$avgComputeSeconds" \
        "$avgWaitSeconds" \
        "$avgFillSeconds" \
        "$avgTotalSeconds" \
        "$checksumValue" \
        "$status" \
        "$logfile" \
        >> "$RESULTS_FILE"
}

runScenario() {
    local algorithm="$1"
    local executable="$2"
    local algorithmArg="$3"
    local scenarioIndex="$4"

    local scenarioName="${SCENARIO_NAMES[$scenarioIndex]}"
    local m="${SCENARIO_M[$scenarioIndex]}"
    local n="${SCENARIO_N[$scenarioIndex]}"
    local k="${SCENARIO_K[$scenarioIndex]}"

    local nodeCount
    local ppn
    local totalRanks
    local nodeList
    local px
    local py
    local output
    local exitStatus
    local metrics
    local checksum
    local avgComputeSeconds
    local avgWaitSeconds
    local avgFillSeconds
    local avgTotalSeconds
    local runsCompleted
    local logfile
    local cannonsDecomp

    for nodeCount in "${NODE_COUNTS[@]}"; do
        echo "Starting ${algorithm} ${scenarioName} on ${nodeCount} node(s)"

        for ppn in "${PPN_VALUES[@]}"; do
            totalRanks=$(( nodeCount * ppn ))
            nodeList="$(buildNodeList "$nodeCount" "$ppn")" || exit 1

            if [[ "$algorithm" == "CANNONS" ]]; then
                cannonsDecomp="$(pickCannonsDecomp "$totalRanks")"

                if [[ -z "$cannonsDecomp" ]]; then
                    read -r px py <<< "$(pickBestDecomp "$totalRanks")"
                    logfile="$(safeLogName "$algorithm" "$scenarioName" "$m" "$nodeCount" "$ppn" "$totalRanks" "$px" "$py")"
                    printf 'SKIP: Cannon\047s requires square Px == Py, but ranks=%d would use %dx%d\n' \
                        "$totalRanks" "$px" "$py" > "$logfile"
                    writeResultRow \
                        "$algorithm" "$scenarioName" "$m" "$n" "$k" \
                        "$nodeCount" "$ppn" "$totalRanks" "$px" "$py" \
                        0 "" "" "" "" "" "SKIP_NON_SQUARE_CANNONS" "$logfile"
                    continue
                fi

                read -r px py <<< "$cannonsDecomp"
            else
                read -r px py <<< "$(pickBestDecomp "$totalRanks")"
            fi

            logfile="$(safeLogName "$algorithm" "$scenarioName" "$m" "$nodeCount" "$ppn" "$totalRanks" "$px" "$py")"

            output="$(
                mpirun.madmpi \
                    -n "$totalRanks" \
                    -nodelist "$nodeList" \
                    "$executable" \
                    "$algorithmArg" \
                    --M "$m" \
                    --N "$n" \
                    --K "$k" \
                    --Px "$px" \
                    --Py "$py" \
                    --runs "$RUNS" \
                    --seed "$SEED" \
                    2>&1
            )"
            exitStatus=$?

            printf '%s\n' "$output" > "$logfile"
            metrics="$(parseMetrics "$output")"
            checksum="$(parseChecksum "$output")"

            if (( exitStatus == 0 )) && [[ -n "$metrics" ]]; then
                read -r avgComputeSeconds avgWaitSeconds avgFillSeconds avgTotalSeconds runsCompleted <<< "$metrics"
                setChecksumCsvValue "$algorithm" "$scenarioName" "$nodeCount" "$ppn" "$checksum"

                writeResultRow \
                    "$algorithm" "$scenarioName" "$m" "$n" "$k" \
                    "$nodeCount" "$ppn" "$totalRanks" "$px" "$py" \
                    "$runsCompleted" "$avgComputeSeconds" "$avgWaitSeconds" \
                    "$avgFillSeconds" "$avgTotalSeconds" "$CHECKSUM_VALUE" "OK" "$logfile"
            else
                if [[ -n "$checksum" ]]; then
                    CHECKSUM_VALUE="$checksum"
                else
                    CHECKSUM_VALUE="NO_CHECKSUM"
                fi

                writeResultRow \
                    "$algorithm" "$scenarioName" "$m" "$n" "$k" \
                    "$nodeCount" "$ppn" "$totalRanks" "$px" "$py" \
                    0 "" "" "" "" "$CHECKSUM_VALUE" "FAIL" "$logfile"
            fi
        done
    done
}

runAlgorithm() {
    local algorithm="$1"
    local executable="$2"
    local algorithmArg="$3"
    local scenarioIndex

    echo "Running all ${algorithm} scenarios"

    for (( scenarioIndex = 0; scenarioIndex < ${#SCENARIO_NAMES[@]}; scenarioIndex++ )); do
        runScenario "$algorithm" "$executable" "$algorithmArg" "$scenarioIndex"
    done
}

runAlgorithm "CANNONS" "$EXECUTABLE" "cannons"

echo "Results written to ${RESULTS_FILE}"
echo "Raw logs written to ${RAW_LOG_DIR}"

