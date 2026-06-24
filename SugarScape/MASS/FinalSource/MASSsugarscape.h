/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: SugarScape                                                      *
 * Library: MASS                                                              *
 *                                                                            *
 * Author: Ahmed Bera Pay                                                     *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef MASS_SUGARSCAPE_H_
#define MASS_SUGARSCAPE_H_

#include "IterationConfig.h"
#include "MASS.h"
#include "SugarChunk.h"

#include <vector>

// Run benchmark from command line arguments
int runSugarScapeBenchmark(
    int argc,
    char* argv[]
);

class MASSsugarscapeBenchmark {
public:
    // Run benchmark from command line arguments
    static int run(
        int argc,
        char* argv[]
    );

private:
    static int gHeight;
    static int gWidth;
    static int gAgentCount;
    static int gTimesteps;
    static int gGrowthRate;
    static int gSugarCapMin;
    static int gSugarCapMax;
    static int gMetabolismMin;
    static int gMetabolismMax;
    static int gVisionMin;
    static int gVisionMax;
    static unsigned int gInitialWealth;

    // Size largest payload sent by exchangeBoundary phases
    static int computeOutMessageBytes(
        int height,
        int width,
        int agentCount,
        int visionMax
    );

    // Bound per-chunk agent report capacity
    static int computeMaxAgentsPerChunk(
        int agentCount,
        int gridWidth,
        int numChunks
    );

    // Gather sugar grid from all MASS places
    static void gatherSugarGrid(
        Places* places,
        int numChunks,
        int height,
        int width,
        std::vector<int>& fullSugar
    );

    // Gather all living agents from all MASS places
    static void gatherAgents(
        Places* places,
        int numChunks,
        int maxAgentsPerChunk,
        std::vector<AgentRec>& allAgents
    );

    // Write optional grid and agent CSV snapshot
    static void writeCsvAt(
        Places* places,
        int numChunks,
        int height,
        int width,
        int stepNumber,
        const char* outputPrefix,
        int maxAgentsPerChunk
    );

    // Build deterministic initial agent payload
    static std::vector<char> buildSeedPayload(
        int agentCount,
        int height,
        int width
    );

    // Build MASS compound iteration config for 1 or more timesteps
    static mass::IterationConfig makeTimestepConfig(
        int iterations,
        bool useAdvanceTimestep,
        const SugarChunkStepArg* stepArg
    );

    // Write run metadata CSV when CSV output is enabled
    static void writeMetadataCsv(
        const char* outputPrefix
    );

    // Write final summary CSV when CSV output is enabled
    static void writeSummaryCsv(
        const char* outputPrefix,
        int liveAgents,
        long long totalWealth,
        int totalSugar
    );

    // Parse command line arguments and reject invalid configuration
    static int parseAndValidateArgs(
        int argc,
        char* argv[],
        int& numProcesses,
        int& numThreads,
        const char*& outputPrefix,
        int& csvInterval,
        bool& debugMode
    );
};

#endif
