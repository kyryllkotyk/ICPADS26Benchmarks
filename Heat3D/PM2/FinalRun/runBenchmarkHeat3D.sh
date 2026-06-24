#!/bin/bash

set -u

RESULTS_FILE="PM2Heat3DResults.csv"
RAW_LOG_DIR="heat3dRawLogs"

DEFAULT_EXECUTABLE="./pm2heat3d"

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

OUTER_RUNS=1
HEAT_RUNS=4

SEED=1
INIT_MIN=0.0
INIT_MAX=100.0
ALPHA=0.4
BETA=0.1

read -r -p "Heat3D executable [${DEFAULT_EXECUTABLE}]: " EXECUTABLE
EXECUTABLE="${EXECUTABLE:-$DEFAULT_EXECUTABLE}"

mkdir -p "$RAW_LOG_DIR"
: > "$RESULTS_FILE"

printf 'scenario,gridX,gridY,gridZ,timesteps,nodes,ppn,ranks,decompX,decompY,decompZ,outerRunsCompleted,heatRunsUsed,avgHaloWallRankMs,avgComputeWallRankMs,avgWallMaxMs,avgHaloMaxMs,avgComputeMaxMs,finalWallRank,finalChecksum,status,logfile\n' \
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
    local scenarioName="$1"
    local gridX="$2"
    local gridY="$3"
    local gridZ="$4"
    local timesteps="$5"
    local nodeCount="$6"
    local ppn="$7"
    local totalRanks="$8"

    printf '%s/%s_X%s_Y%s_Z%s_T%s_nodes%s_ppn%s_ranks%s.log\n' \
        "$RAW_LOG_DIR" \
        "$scenarioName" \
        "$gridX" \
        "$gridY" \
        "$gridZ" \
        "$timesteps" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks"
}

pickBestDecomp() {
    local ranks="$1"
    local bestX=1
    local bestY=1
    local bestZ="$ranks"
    local bestSpread=999999999
    local x
    local y
    local z
    local spread

    for (( x = 1; x <= ranks; x++ )); do
        if (( ranks % x != 0 )); then
            continue
        fi

        for (( y = x; y <= ranks / x; y++ )); do
            if (( (ranks / x) % y != 0 )); then
                continue
            fi

            z=$(( ranks / (x * y) ))

            if (( y > z )); then
                continue
            fi

            spread=$(( z - x ))

            if (( spread < bestSpread )); then
                bestSpread="$spread"
                bestX="$x"
                bestY="$y"
                bestZ="$z"
            fi
        done
    done

    printf '%d %d %d\n' "$bestX" "$bestY" "$bestZ"
}

addFloat() {
    awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", a + b }'
}

avgFloat() {
    awk -v sum="$1" -v count="$2" 'BEGIN { printf "%.3f", sum / count }'
}

parseMetrics() {
    local output="$1"

    printf '%s\n' "$output" | awk '
        /HaloMaxMs=/ && !/^Run 0 / {
            for (i = 1; i <= NF; i++) {
                split($i, pair, "=")

                if (pair[1] == "HaloWallRankMs") {
                    haloWall += pair[2]
                }
                else if (pair[1] == "ComputeWallRankMs") {
                    computeWall += pair[2]
                }
                else if (pair[1] == "WallMaxMs") {
                    wallMax += pair[2]
                }
                else if (pair[1] == "HaloMaxMs") {
                    haloMax += pair[2]
                }
                else if (pair[1] == "ComputeMaxMs") {
                    computeMax += pair[2]
                }
                else if (pair[1] == "WallRank") {
                    wallRank = pair[2]
                }
                else if (pair[1] == "Checksum") {
                    checksum = pair[2]
                }
            }

            n++
        }

        END {
            if (n > 0) {
                printf "%.3f %.3f %.3f %.3f %.3f %s %s %d\n", \
                    haloWall / n, \
                    computeWall / n, \
                    wallMax / n, \
                    haloMax / n, \
                    computeMax / n, \
                    wallRank, \
                    checksum, \
                    n
            }
        }
    '
}

writeResultRow() {
    local scenarioName="$1"
    local gridX="$2"
    local gridY="$3"
    local gridZ="$4"
    local timesteps="$5"
    local nodeCount="$6"
    local ppn="$7"
    local totalRanks="$8"
    local decompX="$9"
    local decompY="${10}"
    local decompZ="${11}"
    local outerRunsCompleted="${12}"
    local heatRunsUsed="${13}"
    local avgHaloWallRankMs="${14}"
    local avgComputeWallRankMs="${15}"
    local avgWallMaxMs="${16}"
    local avgHaloMaxMs="${17}"
    local avgComputeMaxMs="${18}"
    local finalWallRank="${19}"
    local finalChecksum="${20}"
    local status="${21}"
    local logfile="${22}"

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$scenarioName" \
        "$gridX" \
        "$gridY" \
        "$gridZ" \
        "$timesteps" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks" \
        "$decompX" \
        "$decompY" \
        "$decompZ" \
        "$outerRunsCompleted" \
        "$heatRunsUsed" \
        "$avgHaloWallRankMs" \
        "$avgComputeWallRankMs" \
        "$avgWallMaxMs" \
        "$avgHaloMaxMs" \
        "$avgComputeMaxMs" \
        "$finalWallRank" \
        "$finalChecksum" \
        "$status" \
        "$logfile" \
        >> "$RESULTS_FILE"
}

runScenario() {
    local scenarioName="$1"
    local gridX="$2"
    local gridY="$3"
    local gridZ="$4"
    local timesteps="$5"

    local nodeCount
    local ppn
    local totalRanks
    local nodeList
    local logfile
    local runIndex
    local output
    local exitStatus
    local metrics

    local decompX
    local decompY
    local decompZ

    local haloWallMs
    local computeWallMs
    local wallMaxMs
    local haloMaxMs
    local computeMaxMs
    local wallRank
    local checksum
    local heatRunCount

    local outerRunsCompleted
    local heatRunsUsed

    local haloWallSum
    local computeWallSum
    local wallMaxSum
    local haloMaxSum
    local computeMaxSum

    local avgHaloWallRankMs
    local avgComputeWallRankMs
    local avgWallMaxMs
    local avgHaloMaxMs
    local avgComputeMaxMs

    local finalWallRank
    local finalChecksum

    for nodeCount in "${NODE_COUNTS[@]}"; do
        echo "Starting ${scenarioName} on ${nodeCount} node(s)"

        for ppn in "${PPN_VALUES[@]}"; do
            totalRanks=$(( nodeCount * ppn ))
            nodeList="$(buildNodeList "$nodeCount" "$ppn")" || exit 1
            read -r decompX decompY decompZ <<< "$(pickBestDecomp "$totalRanks")"

            logfile="$(safeLogName "$scenarioName" "$gridX" "$gridY" "$gridZ" "$timesteps" "$nodeCount" "$ppn" "$totalRanks")"

            : > "$logfile"

            outerRunsCompleted=0
            heatRunsUsed=0

            haloWallSum=0
            computeWallSum=0
            wallMaxSum=0
            haloMaxSum=0
            computeMaxSum=0

            finalWallRank=""
            finalChecksum=""

            for (( runIndex = 1; runIndex <= OUTER_RUNS; runIndex++ )); do
                {
                    echo "========== RUN ${runIndex}/${OUTER_RUNS} =========="
                    echo "scenario=${scenarioName}"
                    echo "nodes=${nodeCount} ppn=${ppn} ranks=${totalRanks}"
                    echo "decomp=${decompX}x${decompY}x${decompZ}"
                    echo "command: mpirun.madmpi -n ${totalRanks} -nodelist ${nodeList} ${EXECUTABLE} --timesteps ${timesteps} --gridX ${gridX} --gridY ${gridY} --gridZ ${gridZ} --seed ${SEED} --initMin ${INIT_MIN} --initMax ${INIT_MAX} --alpha ${ALPHA} --beta ${BETA} --runs ${HEAT_RUNS} --debug 0 --decompX ${decompX} --decompY ${decompY} --decompZ ${decompZ}"
                    echo
                } >> "$logfile"

                output="$(
                    mpirun.madmpi \
                        -n "$totalRanks" \
                        -nodelist "$nodeList" \
                        "$EXECUTABLE" \
                        --timesteps "$timesteps" \
                        --gridX "$gridX" \
                        --gridY "$gridY" \
                        --gridZ "$gridZ" \
                        --seed "$SEED" \
                        --initMin "$INIT_MIN" \
                        --initMax "$INIT_MAX" \
                        --alpha "$ALPHA" \
                        --beta "$BETA" \
                        --runs "$HEAT_RUNS" \
                        --debug 0 \
                        --decompX "$decompX" \
                        --decompY "$decompY" \
                        --decompZ "$decompZ" \
                        2>&1
                )"

                exitStatus=$?

                printf '%s\n\n' "$output" >> "$logfile"

                metrics="$(parseMetrics "$output")"

                if (( exitStatus == 0 )) && [[ -n "$metrics" ]]; then
                    read -r \
                        haloWallMs \
                        computeWallMs \
                        wallMaxMs \
                        haloMaxMs \
                        computeMaxMs \
                        wallRank \
                        checksum \
                        heatRunCount \
                        <<< "$metrics"

                    haloWallSum="$(addFloat "$haloWallSum" "$haloWallMs")"
                    computeWallSum="$(addFloat "$computeWallSum" "$computeWallMs")"
                    wallMaxSum="$(addFloat "$wallMaxSum" "$wallMaxMs")"
                    haloMaxSum="$(addFloat "$haloMaxSum" "$haloMaxMs")"
                    computeMaxSum="$(addFloat "$computeMaxSum" "$computeMaxMs")"

                    outerRunsCompleted=$(( outerRunsCompleted + 1 ))
                    heatRunsUsed=$(( heatRunsUsed + heatRunCount ))

                    finalWallRank="$wallRank"
                    finalChecksum="$checksum"
                else
                    {
                        echo "RUN ${runIndex} FAILED"
                        echo
                    } >> "$logfile"
                fi
            done

            if (( outerRunsCompleted > 0 )); then
                avgHaloWallRankMs="$(avgFloat "$haloWallSum" "$outerRunsCompleted")"
                avgComputeWallRankMs="$(avgFloat "$computeWallSum" "$outerRunsCompleted")"
                avgWallMaxMs="$(avgFloat "$wallMaxSum" "$outerRunsCompleted")"
                avgHaloMaxMs="$(avgFloat "$haloMaxSum" "$outerRunsCompleted")"
                avgComputeMaxMs="$(avgFloat "$computeMaxSum" "$outerRunsCompleted")"

                if (( outerRunsCompleted == OUTER_RUNS )); then
                    writeResultRow \
                        "$scenarioName" "$gridX" "$gridY" "$gridZ" "$timesteps" \
                        "$nodeCount" "$ppn" "$totalRanks" \
                        "$decompX" "$decompY" "$decompZ" \
                        "$outerRunsCompleted" "$heatRunsUsed" \
                        "$avgHaloWallRankMs" "$avgComputeWallRankMs" "$avgWallMaxMs" "$avgHaloMaxMs" "$avgComputeMaxMs" \
                        "$finalWallRank" "$finalChecksum" \
                        "OK" "$logfile"
                else
                    writeResultRow \
                        "$scenarioName" "$gridX" "$gridY" "$gridZ" "$timesteps" \
                        "$nodeCount" "$ppn" "$totalRanks" \
                        "$decompX" "$decompY" "$decompZ" \
                        "$outerRunsCompleted" "$heatRunsUsed" \
                        "$avgHaloWallRankMs" "$avgComputeWallRankMs" "$avgWallMaxMs" "$avgHaloMaxMs" "$avgComputeMaxMs" \
                        "$finalWallRank" "$finalChecksum" \
                        "PARTIAL_FAIL" "$logfile"
                fi
            else
                writeResultRow \
                    "$scenarioName" "$gridX" "$gridY" "$gridZ" "$timesteps" \
                    "$nodeCount" "$ppn" "$totalRanks" \
                    "$decompX" "$decompY" "$decompZ" \
                    0 0 \
                    "" "" "" "" "" \
                    "" "" \
                    "FAIL" "$logfile"
            fi
        done
    done
}

runScenario "S1Medium" 768 768 768 512
runScenario "S2Large" 832 832 832 512

echo "Results written to ${RESULTS_FILE}"
echo "Raw logs written to ${RAW_LOG_DIR}"
