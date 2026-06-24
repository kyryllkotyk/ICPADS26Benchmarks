#include "MASSsugarscape.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

int MASSsugarscapeBenchmark::gHeight = 0;
int MASSsugarscapeBenchmark::gWidth = 0;
int MASSsugarscapeBenchmark::gAgentCount = 0;
int MASSsugarscapeBenchmark::gTimesteps = 100;
int MASSsugarscapeBenchmark::gGrowthRate = 1;
int MASSsugarscapeBenchmark::gSugarCapMin = 1;
int MASSsugarscapeBenchmark::gSugarCapMax = 4;
int MASSsugarscapeBenchmark::gMetabolismMin = 1;
int MASSsugarscapeBenchmark::gMetabolismMax = 4;
int MASSsugarscapeBenchmark::gVisionMin = 1;
int MASSsugarscapeBenchmark::gVisionMax = 6;
unsigned int MASSsugarscapeBenchmark::gInitialWealth =
    std::numeric_limits<unsigned int>::max();

int MASSsugarscapeBenchmark::computeOutMessageBytes(
    const int height,
    const int width,
    const int agentCount,
    const int visionMax
) {
    (void)height;
    (void)agentCount;

    // Size largest phase payload used by boundary exchange
    const int haloBytes = 2 * visionMax * width * (int)sizeof(int);
    const int candidateSlots = 2 * width;
    const int candidateBytes =
        2 * candidateSlots * (int)sizeof(MoveCandidateRec);
    const int migrationSlots = 2 * width;
    const int migrationBytes =
        2 * migrationSlots * (int)sizeof(AgentRec);

    const int payloadBytes = std::max(
        {
            haloBytes,
            candidateBytes,
            migrationBytes
        }
    );

    return payloadBytes + 4 * (int)sizeof(int) + 64;
}

int MASSsugarscapeBenchmark::computeMaxAgentsPerChunk(
    const int agentCount,
    const int gridWidth,
    const int numChunks
) {
    const int safeChunkCount = std::max(
        numChunks,
        1
    );
    const int averagePerChunk =
        (agentCount + safeChunkCount - 1) / safeChunkCount;

    int maxAgentsPerChunk = std::max(
        1024,
        std::max(
            2 * gridWidth,
            4 * averagePerChunk
        )
    );

    if (agentCount > 0 && maxAgentsPerChunk > 2 * agentCount) {
        maxAgentsPerChunk = 2 * agentCount;
    }

    return maxAgentsPerChunk;
}

void MASSsugarscapeBenchmark::gatherSugarGrid(
    Places* places,
    const int numChunks,
    const int height,
    const int width,
    std::vector<int>& fullSugar
) {
    const int rowsPerChunk = (height + numChunks - 1) / numChunks;
    const int intsPerChunk = 2 + rowsPerChunk * width;
    const int sugarReturnBytes = (int)sizeof(int) * intsPerChunk;

    int* sugarBuffer = (int*)places->callAll(
        "SugarChunk::reportSugar",
        (void**)nullptr,
        0,
        sugarReturnBytes
    );

    fullSugar.assign(
        height * width,
        0
    );

    for (int chunk = 0; chunk < numChunks; chunk++) {
        const int* chunkBuffer = sugarBuffer + chunk * intsPerChunk;
        const int startRow = chunkBuffer[0];
        const int rowCount = chunkBuffer[1];
        const int* data = chunkBuffer + 2;

        for (int row = 0; row < rowCount; row++) {
            for (int col = 0; col < width; col++) {
                fullSugar[(startRow + row) * width + col] =
                    data[row * width + col];
            }
        }
    }

    delete[] sugarBuffer;
}

void MASSsugarscapeBenchmark::gatherAgents(
    Places* places,
    const int numChunks,
    const int maxAgentsPerChunk,
    std::vector<AgentRec>& allAgents
) {
    const int agentReturnBytes =
        (int)sizeof(int) + maxAgentsPerChunk * (int)sizeof(AgentRec);

    char* agentBuffer = (char*)places->callAll(
        "SugarChunk::reportAgents",
        (void**)nullptr,
        0,
        agentReturnBytes
    );

    allAgents.clear();

    for (int chunk = 0; chunk < numChunks; chunk++) {
        const char* chunkBuffer = agentBuffer + chunk * agentReturnBytes;
        const int agentCount = *reinterpret_cast<const int*>(
            chunkBuffer
        );

        if (agentCount < 0) {
            std::cerr
                << "WARNING: reportAgents overflow on chunk "
                << chunk
                << " during CSV export. CSV may be incomplete. "
                << "BENCH_STAT is not affected."
                << std::endl;
            continue;
        }

        const AgentRec* agents = reinterpret_cast<const AgentRec*>(
            chunkBuffer + sizeof(int)
        );

        for (
            int agentIndex = 0;
            agentIndex < agentCount && agentIndex < maxAgentsPerChunk;
            agentIndex++
        ) {
            allAgents.push_back(
                agents[agentIndex]
            );
        }
    }

    delete[] agentBuffer;

    std::sort(
        allAgents.begin(),
        allAgents.end(),
        [](
            const AgentRec& leftAgent,
            const AgentRec& rightAgent
        ) {
            return leftAgent.id < rightAgent.id;
        }
    );
}

void MASSsugarscapeBenchmark::writeCsvAt(
    Places* places,
    const int numChunks,
    const int height,
    const int width,
    const int stepNumber,
    const char* outputPrefix,
    const int maxAgentsPerChunk
) {
    if (outputPrefix == NULL) {
        return;
    }

    std::vector<int> fullSugar;
    std::vector<AgentRec> allAgents;

    gatherSugarGrid(
        places,
        numChunks,
        height,
        width,
        fullSugar
    );
    gatherAgents(
        places,
        numChunks,
        maxAgentsPerChunk,
        allAgents
    );

    std::ostringstream gridPath;
    gridPath
        << outputPrefix
        << "_"
        << std::setfill(
            '0'
        )
        << std::setw(
            6
        )
        << stepNumber
        << "_GRID.csv";

    std::ofstream gridFile(
        gridPath.str()
    );
    gridFile << "row,col,sugarCount\n";

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            gridFile
                << row
                << ","
                << col
                << ","
                << fullSugar[row * width + col]
                << "\n";
        }
    }

    gridFile.close();

    std::ostringstream agentPath;
    agentPath
        << outputPrefix
        << "_"
        << std::setfill(
            '0'
        )
        << std::setw(
            6
        )
        << stepNumber
        << "_AGENT.csv";

    std::ofstream agentFile(
        agentPath.str()
    );
    agentFile << "agentID,row,col,wealth,metabolism,vision,age\n";

    for (const auto& agent : allAgents) {
        agentFile
            << agent.id
            << ","
            << agent.row
            << ","
            << agent.col
            << ","
            << agent.wealth
            << ","
            << agent.metabolism
            << ","
            << agent.vision
            << ","
            << agent.age
            << "\n";
    }

    agentFile.close();
}

std::vector<char> MASSsugarscapeBenchmark::buildSeedPayload(
    const int agentCount,
    const int height,
    const int width
) {
    SugarXoshiro256 metabolismRng(
        2
    );
    SugarXoshiro256 visionRng(
        3
    );
    SugarXoshiro256 coordinateRng(
        4
    );
    SugarXoshiro256 wealthRng(
        100
    );

    std::vector<AgentRec> allAgents(
        agentCount
    );

    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        AgentRec& agent = allAgents[agentIndex];
        agent.id = agentIndex;
        agent.metabolism = (int)metabolismRng.nextInRange(
            gMetabolismMin,
            gMetabolismMax
        );
        agent.vision = (int)visionRng.nextInRange(
            gVisionMin,
            gVisionMax
        );

        if (gInitialWealth == std::numeric_limits<unsigned int>::max()) {
            agent.wealth = (int)wealthRng.nextInRange(
                (uint64_t)(agent.metabolism * 5),
                (uint64_t)(agent.metabolism * 10)
            );
        }
        else {
            agent.wealth = (int)gInitialWealth;
        }

        agent.age = 0;
    }

    std::vector<int> shuffledRows(
        height * width
    );
    std::vector<int> shuffledCols(
        height * width
    );

    int cellIndex = 0;
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            shuffledRows[cellIndex] = row;
            shuffledCols[cellIndex] = col;
            cellIndex++;
        }
    }

    for (int index = (int)shuffledRows.size() - 1; index > 0; index--) {
        const int swapIndex = (int)coordinateRng.nextInRange(
            0,
            (uint64_t)index
        );

        if (swapIndex != index) {
            std::swap(
                shuffledRows[index],
                shuffledRows[swapIndex]
            );
            std::swap(
                shuffledCols[index],
                shuffledCols[swapIndex]
            );
        }
    }

    for (int agentIndex = 0; agentIndex < agentCount; agentIndex++) {
        allAgents[agentIndex].row = shuffledRows[agentIndex];
        allAgents[agentIndex].col = shuffledCols[agentIndex];
    }

    const size_t payloadBytes =
        sizeof(SugarChunkSeedHeader) + agentCount * sizeof(AgentRec);

    std::vector<char> seedBuffer(
        payloadBytes
    );

    SugarChunkSeedHeader* header = reinterpret_cast<SugarChunkSeedHeader*>(
        seedBuffer.data()
    );
    header->numAgents = agentCount;

    std::memcpy(
        seedBuffer.data() + sizeof(SugarChunkSeedHeader),
        allAgents.data(),
        agentCount * sizeof(AgentRec)
    );

    return seedBuffer;
}

mass::IterationConfig MASSsugarscapeBenchmark::makeTimestepConfig(
    const int iterations,
    const bool useAdvanceTimestep,
    const SugarChunkStepArg* stepArg
) {
    mass::IterationConfig config;
    config.iterations(
        iterations
    );

    if (useAdvanceTimestep) {
        config.placeCompute(
            "SugarChunk::advanceTimestep"
        );
    }
    else {
        config.placeCompute(
            "SugarChunk::setTimestep",
            const_cast<void*>(
                static_cast<const void*>(
                    stepArg
                )
            ),
            sizeof(SugarChunkStepArg)
        );
    }

    config
        .placeCompute(
            "SugarChunk::growSugar"
        )
        .placeCompute(
            "SugarChunk::packHaloSugar"
        )
        .placeExchangeBoundary()
        .placeCompute(
            "SugarChunk::computeTargets"
        )
        .placeCompute(
            "SugarChunk::packMoveCandidates"
        )
        .placeExchangeBoundary()
        .placeCompute(
            "SugarChunk::resolveConflictsLocal"
        )
        .placeCompute(
            "SugarChunk::packWinnerReplies"
        )
        .placeExchangeBoundary()
        .placeCompute(
            "SugarChunk::applyMovesAndTakeSugar"
        )
        .placeCompute(
            "SugarChunk::packMigrations"
        )
        .placeExchangeBoundary()
        .placeCompute(
            "SugarChunk::integrateMigrations"
        )
        .placeCompute(
            "SugarChunk::sortLocalAgents"
        )
        .placeCompute(
            "SugarChunk::takeSugarAtCurrent"
        )
        .placeCompute(
            "SugarChunk::metabolismAndDeath"
        );

    return config;
}

void MASSsugarscapeBenchmark::writeMetadataCsv(
    const char* outputPrefix
) {
    if (outputPrefix == NULL) {
        return;
    }

    std::ostringstream metadataPath;
    metadataPath << outputPrefix << "_META.csv";

    std::ofstream metadataFile(
        metadataPath.str()
    );

    metadataFile
        << "Height,width,timesteps,growthRate,sugarCapacityMin,"
        << "sugarCapacityMax,initialAgentCount,metabolismMin,"
        << "metabolismMax,visionMin,visionMax,initialWealth,"
        << "seedSugar,seedMetabolism,seedVision,seedCoord,seedWealth\n";

    metadataFile
        << gHeight
        << ","
        << gWidth
        << ","
        << gTimesteps
        << ","
        << gGrowthRate
        << ","
        << gSugarCapMin
        << ","
        << gSugarCapMax
        << ","
        << gAgentCount
        << ","
        << gMetabolismMin
        << ","
        << gMetabolismMax
        << ","
        << gVisionMin
        << ","
        << gVisionMax
        << ",";

    if (gInitialWealth == std::numeric_limits<unsigned int>::max()) {
        metadataFile << "UINT_MAX";
    }
    else {
        metadataFile << gInitialWealth;
    }

    metadataFile << ",1,2,3,4,100\n";
    metadataFile.close();
}

void MASSsugarscapeBenchmark::writeSummaryCsv(
    const char* outputPrefix,
    const int liveAgents,
    const long long totalWealth,
    const int totalSugar
) {
    if (outputPrefix == NULL) {
        return;
    }

    std::ostringstream summaryPath;
    summaryPath
        << outputPrefix
        << "_"
        << std::setfill(
            '0'
        )
        << std::setw(
            6
        )
        << gTimesteps
        << "_TIMESTEPCONSTANTS.csv";

    std::ofstream summaryFile(
        summaryPath.str()
    );

    summaryFile
        << "timestep,liveAgentCount,totalWealth,"
        << "averageWealth,totalSugarOnGrid\n";

    summaryFile
        << gTimesteps
        << ","
        << liveAgents
        << ","
        << totalWealth
        << ","
        << (
            liveAgents > 0
                ? (double)totalWealth / liveAgents
                : 0
        )
        << ","
        << totalSugar
        << "\n";

    summaryFile.close();
}

int MASSsugarscapeBenchmark::parseAndValidateArgs(
    const int argc,
    char* argv[],
    int& numProcesses,
    int& numThreads,
    const char*& outputPrefix,
    int& csvInterval,
    bool& debugMode
) {
    if (argc < 10) {
        std::cerr
            << "usage: ./main username password machinefile port nProc nThr "
            << "height width agentCount [timesteps [outputPrefix "
            << "[csvInterval [growthRate [metabolismMax [--debug]]]]]]]"
            << std::endl;
        return -1;
    }

    numProcesses = std::atoi(
        argv[5]
    );
    numThreads = std::atoi(
        argv[6]
    );
    gHeight = std::atoi(
        argv[7]
    );
    gWidth = std::atoi(
        argv[8]
    );
    gAgentCount = std::atoi(
        argv[9]
    );
    gTimesteps = (
        argc > 10
            ? std::atoi(
                argv[10]
            )
            : 100
    );
    outputPrefix = (
        argc > 11
            ? argv[11]
            : NULL
    );
    csvInterval = (
        argc > 12
            ? std::atoi(
                argv[12]
            )
            : 1
    );

    if (csvInterval <= 0) {
        csvInterval = 1;
    }

    if (argc > 13) {
        gGrowthRate = std::atoi(
            argv[13]
        );
    }

    if (argc > 14) {
        gMetabolismMax = std::atoi(
            argv[14]
        );
    }

    debugMode = false;
    for (int argIndex = 15; argIndex < argc; argIndex++) {
        if (std::strcmp(
            argv[argIndex],
            "--debug"
        ) == 0) {
            debugMode = true;
        }
    }

    const long long cellCount =
        (long long)gHeight * (long long)gWidth;

    if (
        numProcesses <= 0 ||
        numThreads <= 0 ||
        gHeight <= 0 ||
        gWidth <= 0 ||
        gAgentCount < 0 ||
        gTimesteps < 0 ||
        gAgentCount > cellCount
    ) {
        std::cerr << "Invalid MASS or SugarScape configuration" << std::endl;
        return -1;
    }

    return 0;
}

int MASSsugarscapeBenchmark::run(
    int argc,
    char* argv[]
) {
    int numProcesses = 0;
    int numThreads = 0;
    const char* outputPrefix = NULL;
    int csvInterval = 1;
    bool debugMode = false;

    if (parseAndValidateArgs(
        argc,
        argv,
        numProcesses,
        numThreads,
        outputPrefix,
        csvInterval,
        debugMode
    ) != 0) {
        return -1;
    }

    const int numChunks = std::max(
        numProcesses * numThreads,
        1
    );

    if (numChunks > 1 && gHeight / numChunks < 2 * gVisionMax) {
        std::cerr
            << "Error: partition height "
            << gHeight / numChunks
            << " < 2*visionMax="
            << 2 * gVisionMax
            << ". Reduce nProc*nThr or increase gridHeight."
            << std::endl;
        return -1;
    }

    auto totalStartTime = std::chrono::steady_clock::now();

    char* massArgs[4] = {
        argv[1],
        argv[2],
        argv[3],
        argv[4]
    };

    MASS::init(
        massArgs,
        numProcesses,
        numThreads
    );

    const int placesHandle = 1;
    const int outMessageBytes = computeOutMessageBytes(
        gHeight,
        gWidth,
        gAgentCount,
        gVisionMax
    );
    const int maxAgentsPerChunk = computeMaxAgentsPerChunk(
        gAgentCount,
        gWidth,
        numChunks
    );
    const int rowsPerChunk = (gHeight + numChunks - 1) / numChunks;
    const int reportSugarInts = 2 + rowsPerChunk * gWidth;
    const int reportAgentsBytes =
        (int)sizeof(int) + maxAgentsPerChunk * (int)sizeof(AgentRec);

    SugarChunkConfig config{};
    config.height = gHeight;
    config.width = gWidth;
    config.numChunks = numChunks;
    config.agentCount = gAgentCount;
    config.growthRate = gGrowthRate;
    config.sugarCapMin = gSugarCapMin;
    config.sugarCapMax = gSugarCapMax;
    config.metabolismMin = gMetabolismMin;
    config.metabolismMax = gMetabolismMax;
    config.visionMin = gVisionMin;
    config.visionMax = gVisionMax;
    config.initialWealth = gInitialWealth;
    config.outMessageBytes = outMessageBytes;
    config.reportAgentsBytes = reportAgentsBytes;
    config.reportSugarInts = reportSugarInts;

    Places* chunks = new Places(
        placesHandle,
        "SugarChunk",
        1,
        &config,
        sizeof(config),
        1,
        numChunks
    );

    chunks->callAll(
        "SugarChunk::init",
        &config,
        sizeof(config)
    );

    {
        std::vector<char> seedBuffer = buildSeedPayload(
            gAgentCount,
            gHeight,
            gWidth
        );

        chunks->callAll(
            "SugarChunk::seedInitialAgents",
            seedBuffer.data(),
            seedBuffer.size()
        );

        std::vector<char>().swap(
            seedBuffer
        );
    }

    chunks->callAll(
        "SugarChunk::sortLocalAgents"
    );

    writeMetadataCsv(
        outputPrefix
    );

    if (outputPrefix != NULL) {
        writeCsvAt(
            chunks,
            numChunks,
            gHeight,
            gWidth,
            0,
            outputPrefix,
            maxAgentsPerChunk
        );
    }

    auto simulationStartTime = std::chrono::steady_clock::now();

    const bool perStepDriver = debugMode || outputPrefix != NULL;

    if (perStepDriver) {
        for (int step = 0; step < gTimesteps; step++) {
            SugarChunkStepArg stepArg{
                step
            };

            MASS::runIterations(
                placesHandle,
                makeTimestepConfig(
                    1,
                    false,
                    &stepArg
                )
            );

            if (debugMode) {
                SugarChunkStats* stats = (SugarChunkStats*)chunks->callAll(
                    "SugarChunk::reportStats",
                    (void**)nullptr,
                    0,
                    sizeof(SugarChunkStats)
                );

                int liveAgents = 0;
                long long totalWealth = 0;
                int totalSugar = 0;

                for (int chunk = 0; chunk < numChunks; chunk++) {
                    liveAgents += stats[chunk].localAgentCount;
                    totalWealth += stats[chunk].localWealth;
                    totalSugar += stats[chunk].localSugar;
                }

                delete[] stats;

                std::cerr
                    << "DEBUG step="
                    << step + 1
                    << " agents="
                    << liveAgents
                    << " wealth="
                    << totalWealth
                    << " sugar="
                    << totalSugar
                    << std::endl;
            }

            const int stepNumber = step + 1;
            if (
                outputPrefix != NULL &&
                (
                    stepNumber == gTimesteps ||
                    stepNumber % csvInterval == 0
                )
            ) {
                writeCsvAt(
                    chunks,
                    numChunks,
                    gHeight,
                    gWidth,
                    stepNumber,
                    outputPrefix,
                    maxAgentsPerChunk
                );
            }
        }
    }
    else {
        MASS::runIterations(
            placesHandle,
            makeTimestepConfig(
                gTimesteps,
                true,
                NULL
            )
        );
    }

    auto totalEndTime = std::chrono::steady_clock::now();

    SugarChunkStats* finalStats = (SugarChunkStats*)chunks->callAll(
        "SugarChunk::reportStats",
        (void**)nullptr,
        0,
        sizeof(SugarChunkStats)
    );

    int liveAgents = 0;
    long long totalWealth = 0;
    int totalSugar = 0;

    for (int chunk = 0; chunk < numChunks; chunk++) {
        liveAgents += finalStats[chunk].localAgentCount;
        totalWealth += finalStats[chunk].localWealth;
        totalSugar += finalStats[chunk].localSugar;
    }

    delete[] finalStats;

    std::cout
        << "BENCH_STAT agents="
        << liveAgents
        << " wealth="
        << totalWealth
        << " sugar="
        << totalSugar
        << std::endl;

    writeSummaryCsv(
        outputPrefix,
        liveAgents,
        totalWealth,
        totalSugar
    );

    const long initializationMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            simulationStartTime - totalStartTime
        ).count();
    const long simulationMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            totalEndTime - simulationStartTime
        ).count();
    const long totalMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            totalEndTime - totalStartTime
        ).count();

    std::cout
        << "INIT_MS="
        << initializationMs
        << std::endl;
    std::cout
        << "WALL_TIME_MS="
        << simulationMs
        << std::endl;
    std::cout
        << "TOTAL_MS="
        << totalMs
        << std::endl;

    MASS::finish();
    return 0;
}

int runSugarScapeBenchmark(
    int argc,
    char* argv[]
) {
    return MASSsugarscapeBenchmark::run(
        argc,
        argv
    );
}
