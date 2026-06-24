#!/bin/bash

set -u

RESULTS_FILE="PM2SugarScapeResults.csv"
RAW_LOG_DIR="sugarscapeRawLogs"

DEFAULT_EXECUTABLE="./pm2sugarscape"

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

SCENARIO_NAME="SugarScape10k10m"
GRID_HEIGHT=10000
GRID_WIDTH=10000
AGENT_COUNT=10000000
TIMESTEPS=100

read -r -p "SugarScape executable [${DEFAULT_EXECUTABLE}]: " EXECUTABLE
EXECUTABLE="${EXECUTABLE:-$DEFAULT_EXECUTABLE}"

mkdir -p "$RAW_LOG_DIR"
: > "$RESULTS_FILE"

printf 'scenario,height,width,agents,timesteps,nodes,ppn,ranks,runsCompleted,avgInitMs,avgWallTimeMs,avgTotalMs,finalLiveAgents,finalWealth,finalSugar,status,logfile\n' \
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

safeLogName() {
    local nodeCount="$1"
    local ppn="$2"
    local totalRanks="$3"

    printf '%s/%s_H%s_W%s_A%s_T%s_nodes%s_ppn%s_ranks%s.log\n' \
        "$RAW_LOG_DIR" \
        "$SCENARIO_NAME" \
        "$GRID_HEIGHT" \
        "$GRID_WIDTH" \
        "$AGENT_COUNT" \
        "$TIMESTEPS" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks"
}

parseMetrics() {
    local output="$1"

    printf '%s\n' "$output" | awk '
        /^BENCH_STAT / {
            for (i = 1; i <= NF; i++) {
                split($i, pair, "=")

                if (pair[1] == "agents") {
                    agents = pair[2]
                }
                else if (pair[1] == "wealth") {
                    wealth = pair[2]
                }
                else if (pair[1] == "sugar") {
                    sugar = pair[2]
                }
            }
        }

        /^INIT_MS=/ {
            split($0, pair, "=")
            initMs = pair[2]
        }

        /^WALL_TIME_MS=/ {
            split($0, pair, "=")
            wallMs = pair[2]
        }

        /^TOTAL_MS=/ {
            split($0, pair, "=")
            totalMs = pair[2]
        }

        END {
            if (initMs != "" && wallMs != "" && totalMs != "" && agents != "" && wealth != "" && sugar != "") {
                printf "%s %s %s %s %s %s\n", initMs, wallMs, totalMs, agents, wealth, sugar
            }
        }
    '
}

writeResultRow() {
    local nodeCount="$1"
    local ppn="$2"
    local totalRanks="$3"
    local runsCompleted="$4"
    local avgInitMs="$5"
    local avgWallTimeMs="$6"
    local avgTotalMs="$7"
    local finalLiveAgents="$8"
    local finalWealth="$9"
    local finalSugar="${10}"
    local status="${11}"
    local logfile="${12}"

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$SCENARIO_NAME" \
        "$GRID_HEIGHT" \
        "$GRID_WIDTH" \
        "$AGENT_COUNT" \
        "$TIMESTEPS" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks" \
        "$runsCompleted" \
        "$avgInitMs" \
        "$avgWallTimeMs" \
        "$avgTotalMs" \
        "$finalLiveAgents" \
        "$finalWealth" \
        "$finalSugar" \
        "$status" \
        "$logfile" \
        >> "$RESULTS_FILE"
}

runScenario() {
    local nodeCount
    local ppn
    local totalRanks
    local nodeList
    local logfile
    local runIndex
    local output
    local exitStatus
    local metrics
    local initMs
    local wallMs
    local totalMs
    local liveAgents
    local wealth
    local sugar
    local runsCompleted
    local initSum
    local wallSum
    local totalSum
    local avgInitMs
    local avgWallTimeMs
    local avgTotalMs

    for nodeCount in "${NODE_COUNTS[@]}"; do
        echo "Starting ${SCENARIO_NAME} on ${nodeCount} node(s)"

        for ppn in "${PPN_VALUES[@]}"; do
            totalRanks=$(( nodeCount * ppn ))
            nodeList="$(buildNodeList "$nodeCount" "$ppn")" || exit 1
            logfile="$(safeLogName "$nodeCount" "$ppn" "$totalRanks")"

            : > "$logfile"

            runsCompleted=0
            initSum=0
            wallSum=0
            totalSum=0
            liveAgents=""
            wealth=""
            sugar=""

            for (( runIndex = 1; runIndex <= RUNS; runIndex++ )); do
                {
                    echo "========== RUN ${runIndex}/${RUNS} =========="
                    echo "nodes=${nodeCount} ppn=${ppn} ranks=${totalRanks}"
                    echo "command: mpirun.madmpi -n ${totalRanks} -nodelist ${nodeList} ${EXECUTABLE} ${GRID_HEIGHT} ${GRID_WIDTH} ${AGENT_COUNT} ${TIMESTEPS}"
                    echo
                } >> "$logfile"

                output="$(
                    mpirun.madmpi \
                        -n "$totalRanks" \
                        -nodelist "$nodeList" \
                        "$EXECUTABLE" \
                        "$GRID_HEIGHT" \
                        "$GRID_WIDTH" \
                        "$AGENT_COUNT" \
                        "$TIMESTEPS" \
                        2>&1
                )"

                exitStatus=$?

                printf '%s\n\n' "$output" >> "$logfile"

                metrics="$(parseMetrics "$output")"

                if (( exitStatus == 0 )) && [[ -n "$metrics" ]]; then
                    read -r initMs wallMs totalMs liveAgents wealth sugar <<< "$metrics"

                    initSum=$(( initSum + initMs ))
                    wallSum=$(( wallSum + wallMs ))
                    totalSum=$(( totalSum + totalMs ))
                    runsCompleted=$(( runsCompleted + 1 ))
                else
                    {
                        echo "RUN ${runIndex} FAILED"
                        echo
                    } >> "$logfile"
                fi
            done

            if (( runsCompleted > 0 )); then
                avgInitMs="$(awk -v sum="$initSum" -v count="$runsCompleted" 'BEGIN { printf "%.3f", sum / count }')"
                avgWallTimeMs="$(awk -v sum="$wallSum" -v count="$runsCompleted" 'BEGIN { printf "%.3f", sum / count }')"
                avgTotalMs="$(awk -v sum="$totalSum" -v count="$runsCompleted" 'BEGIN { printf "%.3f", sum / count }')"

                if (( runsCompleted == RUNS )); then
                    writeResultRow \
                        "$nodeCount" "$ppn" "$totalRanks" \
                        "$runsCompleted" "$avgInitMs" "$avgWallTimeMs" "$avgTotalMs" \
                        "$liveAgents" "$wealth" "$sugar" \
                        "OK" "$logfile"
                else
                    writeResultRow \
                        "$nodeCount" "$ppn" "$totalRanks" \
                        "$runsCompleted" "$avgInitMs" "$avgWallTimeMs" "$avgTotalMs" \
                        "$liveAgents" "$wealth" "$sugar" \
                        "PARTIAL_FAIL" "$logfile"
                fi
            else
                writeResultRow \
                    "$nodeCount" "$ppn" "$totalRanks" \
                    0 "" "" "" "" "" "" \
                    "FAIL" "$logfile"
            fi
        done
    done
}

runScenario

echo "Results written to ${RESULTS_FILE}"
echo "Raw logs written to ${RAW_LOG_DIR}"
