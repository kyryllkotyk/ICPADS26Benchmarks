#!/bin/bash

set -u

RESULTS_FILE="PM2MotifSearchResults.csv"
RAW_LOG_DIR="motifSearchRawLogsEdgefile"
TMP_CSV_DIR="motifSearchTmpCsvEdgefile"

DEFAULT_EXECUTABLE="./pm2motifsearch"
DEFAULT_EDGE_FILE="triangle250edgelist.edgelist"

MOTIF_NAME="triangle"
VERTICES="0"
EDGE_PROB="0.0011"
SEED="42"
DEBUG="0"
RUNS="3"

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

read -r -p "MotifSearch executable [${DEFAULT_EXECUTABLE}]: " MOTIF_EXECUTABLE
MOTIF_EXECUTABLE="${MOTIF_EXECUTABLE:-$DEFAULT_EXECUTABLE}"

read -r -p "Edge file [${DEFAULT_EDGE_FILE}]: " EDGE_FILE
EDGE_FILE="${EDGE_FILE:-$DEFAULT_EDGE_FILE}"

if [[ ! -x "$MOTIF_EXECUTABLE" ]]; then
    echo "ERROR: executable not found or not executable: $MOTIF_EXECUTABLE" >&2
    exit 1
fi

if [[ ! -f "$EDGE_FILE" ]]; then
    echo "ERROR: edge file not found: $EDGE_FILE" >&2
    exit 1
fi

mkdir -p "$RAW_LOG_DIR"
mkdir -p "$TMP_CSV_DIR"

: > "$RESULTS_FILE"

printf 'implementation,motif,vertices,edge_prob,seed,total_processes,nodes_used,processes_per_node,edges,raw_embeddings,automorphisms,instance_count,search_time_ms,total_wall_time_ms,status,logfile\n' \
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

    printf '%s/motif_edgefile_nodes%s_ppn%s_ranks%s.log\n' \
        "$RAW_LOG_DIR" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks"
}

safeTmpCsvName() {
    local nodeCount="$1"
    local ppn="$2"
    local totalRanks="$3"

    printf '%s/motif_edgefile_nodes%s_ppn%s_ranks%s.csv\n' \
        "$TMP_CSV_DIR" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks"
}

safeBatchCsvName() {
    local nodeCount="$1"
    local ppn="$2"
    local totalRanks="$3"

    printf '%s/motif_edgefile_nodes%s_ppn%s_ranks%s_batch.csv\n' \
        "$TMP_CSV_DIR" \
        "$nodeCount" \
        "$ppn" \
        "$totalRanks"
}

writeBatchCsv() {
    local batchCsv="$1"

    {
        printf 'motif,vertices,edge_prob,seed,input,debug,runs\n'
        printf '%s,%s,%s,%s,%s,%s,%s\n' \
            "$MOTIF_NAME" \
            "$VERTICES" \
            "$EDGE_PROB" \
            "$SEED" \
            "$EDGE_FILE" \
            "$DEBUG" \
            "$RUNS"
    } > "$batchCsv"
}

appendTmpCsvRows() {
    local tmpCsv="$1"
    local logfile="$2"

    if [[ ! -s "$tmpCsv" ]]; then
        return 1
    fi

    awk \
        -F',' \
        -v OFS=',' \
        -v logfile="$logfile" '
        NR == 1 {
            next
        }

        NF > 0 {
            $NF = logfile
            print
        }
    ' "$tmpCsv" >> "$RESULTS_FILE"

    return 0
}

writeFailRow() {
    local totalRanks="$1"
    local nodeCount="$2"
    local ppn="$3"
    local logfile="$4"

    printf 'pm2,%s,%s,%s,%s,%s,%s,%s,,,,,,,%s,%s\n' \
        "$MOTIF_NAME" \
        "$VERTICES" \
        "$EDGE_PROB" \
        "$SEED" \
        "$totalRanks" \
        "$nodeCount" \
        "$ppn" \
        "FAIL" \
        "$logfile" \
        >> "$RESULTS_FILE"
}

runConfig() {
    local nodeCount="$1"
    local ppn="$2"

    local totalRanks=$(( nodeCount * ppn ))
    local nodeList
    local logfile
    local tmpCsv
    local batchCsv
    local exitStatus

    nodeList="$(buildNodeList "$nodeCount" "$ppn")" || exit 1
    logfile="$(safeLogName "$nodeCount" "$ppn" "$totalRanks")"
    tmpCsv="$(safeTmpCsvName "$nodeCount" "$ppn" "$totalRanks")"
    batchCsv="$(safeBatchCsvName "$nodeCount" "$ppn" "$totalRanks")"

    rm -f "$tmpCsv"
    rm -f "$batchCsv"

    writeBatchCsv "$batchCsv"

    echo "Running MotifSearch edgefile nodes=${nodeCount} ppn=${ppn} ranks=${totalRanks}"

    mpirun.madmpi \
        -n "$totalRanks" \
        -nodelist "$nodeList" \
        "$MOTIF_EXECUTABLE" \
        --batch-file "$batchCsv" \
        --output "$tmpCsv" \
        --nodes-used "$nodeCount" \
        --processes-per-node "$ppn" \
        > "$logfile" 2>&1

    exitStatus=$?

    if (( exitStatus == 0 )) && appendTmpCsvRows "$tmpCsv" "$logfile"; then
        echo "OK nodes=${nodeCount} ppn=${ppn} ranks=${totalRanks}"
    else
        echo "FAIL nodes=${nodeCount} ppn=${ppn} ranks=${totalRanks}; see $logfile"
        writeFailRow "$totalRanks" "$nodeCount" "$ppn" "$logfile"
    fi
}

echo "Starting MotifSearch edgefile benchmark run"
echo "Executable: $MOTIF_EXECUTABLE"
echo "Edge file:  $EDGE_FILE"
echo "Motif:      $MOTIF_NAME"
echo "Runs:       $RUNS"
echo "Results:    $RESULTS_FILE"
echo "Raw logs:   $RAW_LOG_DIR"

for nodeCount in "${NODE_COUNTS[@]}"; do
    echo "=== ${nodeCount} node(s) ==="

    for ppn in "${PPN_VALUES[@]}"; do
        runConfig "$nodeCount" "$ppn"
    done
done

echo "Results written to ${RESULTS_FILE}"
echo "Raw logs written to ${RAW_LOG_DIR}"
echo "Temporary CSVs written to ${TMP_CSV_DIR}"
