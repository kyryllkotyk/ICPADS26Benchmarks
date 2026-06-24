/******************************************************************************
 *                                                                            *
 * ICPADS26 Benchmarks Collection                                             *
 *                                                                            *
 * Benchmark: Bail-In Bail-Out                                                *
 * Library: PM2/MPI                                                           *
 *                                                                            *
 * Author: Kyryll Kotyk                                                       *
 * Faculty Advisor: Munehiro Fukuda                                           *
 * Code Finalization: Kyryll Kotyk                                            *
 *                                                                            *
 *****************************************************************************/

#ifndef PM2_BAIL_IN_BAIL_OUT_
#define PM2_BAIL_IN_BAIL_OUT_

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

using namespace std;

// RNG stream IDs for firm production and shock parameters
#define STREAM_FIRM_SHOCK_PERCENT 1000ULL
#define STREAM_FIRM_PROFIT_PERCENT 1001ULL
#define STREAM_FIRM_OPERATING_COST_PERCENT 1002ULL
#define STREAM_FIRM_SHOCK_PROFIT_MULTIPLIER_PERCENT 1003ULL
#define STREAM_FIRM_INITIAL_LIQUIDITY_MULTIPLIER_PERCENT 1004ULL
#define STREAM_FIRM_RARE_EVENT_PROBABILITY 1005ULL
#define STREAM_FIRM_RARE_EVENT_IMPACT_PERCENT 1006ULL
#define STREAM_FIRM_RARE_EVENT_SIGN 1007ULL
#define STREAM_SHOCK_SIGN 1008ULL

#define STREAM_BANK_GRAPH_BUILD 1ULL
#define STREAM_WORKER_BANK_ASSIGN 2ULL
#define STREAM_WORKER_FIRM_ASSIGN 3ULL
#define STREAM_INTEREST_RATE_SELECTION 4ULL
#define STREAM_SHOCK_GENERATION 5ULL
#define STREAM_PROFIT_MULTIPLIER_GENERATION 6ULL
#define STREAM_FIRM_BANK_GRAPH_BUILD 7ULL
#define STREAM_INTERBANK_LENDER_SAMPLING 8ULL

class BailInBailOut
{
public:
    /**
    * @brief Executes the Bail-In Bail-Out financial contagion benchmark
    *
    * @details This function runs the complete simulation of a financial system
    * containing banks, firms, and workers. The simulation is distributed across
    * MPI ranks, where each rank owns a part of the workers, firms, and banks.
    * The benchmark models firm production, worker wage deposits,
    * bank lending, interbank loans, and financial distress policies
    *
    * For each run, two policy modes are executed:
    *
    * Bail-In (interbank debts are forgiven by the lending banks)
    * Bail-Out (external liquidity is injected into insolvent banks)
    *
    * Each run consists of multiple timesteps. During every timestep the following
    * phases are executed:
    *
    * A: Firm production, shocks, profit realization, and loan repayment
    * B: Worker consumption and wage deposits
    * C: Worker deposits and firm repayments transferred to banks
    * D: Firms request loans from banks and banks process requests
    * E: Banks repay interbank loans, accrue interest, and borrow if insolvent
    * F: Insolvency detection and application of bail-in or bail-out policy
    *
    * The benchmark is designed to stress distributed systems through:
    *
    * Irregular message passing (loan requests and acceptances)
    * Graph-based communication (interbank graph network)
    * Distributed state updates
    * Deterministic randomization
    *
    * All randomness is deterministically seeded so that identical results are
    * produced regardless of MPI system configuration.
    *
    * @param runs Number of independent simulation runs to perform
    * @param timesteps Number of timesteps within each run
    * @param baseSeed Base seed used for random number generation
    * The seed is mixed with run number, timestep, IDs, and stream tags for
    * reproducibility
    *
    * @param bankCountTotal Total number of banks in the simulation
    * @param firmCountTotal Total number of firms in the simulation
    * @param workerCountTotal Total number of workers in the simulation
    * @param bankWorkerCount Number of employees assigned to each bank.
    * Determines wage payments made by banks every timestep
    *
    * @param initialBankLiquidity Initial liquidity assigned to every bank
    * @param initialFirmLiquidity Initial liquidity assigned to every firm
    * @param initialProductionCost Initial production cost assigned to firms.
    * Production cost changes over time due to shocks and profit multipliers
    * @param wage Wage paid to each firm worker per timestep
    * @param bankEmployeeWage Wage paid to each bank employee per timestep
    *
    * @param wageConsumptionPercent Percentage of wages consumed by workers.
    * The remaining portion is deposited into worker bank accounts
    * @param profitMultiplierMin Minimum profit multiplier applied to firm
    * production output (in percent)
    * @param profitMultiplierMax Maximum profit multiplier applied to firm
    * production output (in percent)
    * @param shockMultiplierMin Minimum shock multiplier applied to firm
    * production cost (in percent)
    * @param shockMultiplierMax Maximum shock multiplier applied to firm
    * production cost (in percent)
    * @param minInterestRate Minimum interest rate applied to loans (in percent)
    * @param maxInterestRate Maximum interest rate applied to loans (in percent)
    *
    * @param interbankDensity Density of the interbank lending graph (percentage
    * of possible edges used when constructing neighbor lists)
    * @param maxInterbankLenderSamplingK Maximum number of lenders sampled when a
    * bank tries borrowing from neighbors in the interbank network
    * @param maxInterbankLoanPercent Maximum percentage of a bank's liquidity that
    * can be loaned to another bank during interbank lending
    * @param maxFirmLoanPercent Maximum percentage of a bank's liquidity that may
    * be lent to firms in a single timestep
    * @param firmLenderDegree Number of banks each firm is
    * connected to in the firm-bank graph
    * @param firmRepayPercent Percentage of outstanding firm debt repaid each
    * timestep when liquidity permits
    * @param bankRepayPercent Percentage of outstanding interbank debt repaid
    * each timestep when liquidity permits
    *
    * @param interventionDelay How many timesteps to wait before applying
    * Bail-In or Bail-Out policy to the insolvent bank
    *
    * @param bailInCoveragePercent Percentage of a bank's debts that
    * will be forgiven during a bail-in intervention
    *
    * @param bailOutCoveragePercent Percentage of a bank's deficit covered by
    * external liquidity injection during a bail-out intervention
    *
    */
    void runBenchmarkInAndOut(
        /* Simulation Specifiers */
        const unsigned int runs,
        const unsigned int timesteps,
        const uint64_t baseSeed,

        /* Structural Parameters */
        const unsigned int bankCountTotal,
        const unsigned int firmCountTotal,
        const unsigned int workerCountTotal,
        // How many workers each bank has
        const unsigned int bankWorkerCount,

        /* System Parameters */
        const unsigned int initialBankLiquidity,
        const unsigned int initialFirmLiquidity,
        const unsigned int initialProductionCost,
        const unsigned int wage,
        const unsigned int bankEmployeeWage,

        /* Multipliers (in %, 95 = 0.95 multiplier) */
        unsigned short wageConsumptionPercent,
        unsigned short profitMultiplierMin,
        unsigned short profitMultiplierMax,
        unsigned short shockMultiplierMin,
        unsigned short shockMultiplierMax,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        unsigned short firmShockMinPercent,
        unsigned short firmShockMaxPercent,
        unsigned short firmProfitMinPercent,
        unsigned short firmProfitMaxPercent,
        unsigned short firmOperatingCostMinPercent,
        unsigned short firmOperatingCostMaxPercent,
        unsigned short firmShockProfitMultiplierMinPercent,
        unsigned short firmShockProfitMultiplierMaxPercent,
        unsigned short firmInitialLiquidityMultiplierMinPercent,
        unsigned short firmInitialLiquidityMultiplierMaxPercent,
        unsigned short firmRareEventProbabilityPercent,
        unsigned short firmRareEventImpactMinPercent,
        unsigned short firmRareEventImpactMaxPercent,

        /* Banking Interaction Parameters */
        unsigned short interbankDensity,
        const unsigned short maxInterbankLenderSamplingK,
        // Maximum percentage of the bank's liquidity it can offer as a loan
        unsigned short maxInterbankLoanPercent,
        unsigned short maxFirmLoanPercent,
        unsigned short firmLenderDegree,
        unsigned short firmRepayPercent,
        unsigned short bankRepayPercent,

        unsigned int interventionDelay,

        /* Bail-In Parameters */
        unsigned short bailInCoveragePercent,

        /* Bail-Out Parameters */
        unsigned short bailOutCoveragePercent
    );

private:
    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_rng.cpp                                               *
    **************************************************************************/

    //===============================STRUCTURES================================

    struct Xoshiro256 {
        uint64_t s[4];

        explicit Xoshiro256(
            uint64_t seed
        );

        static uint64_t splitmix64(
            uint64_t& x
        );
        static uint64_t rotl(
            uint64_t x,
            int k
        );

        uint64_t next();
        uint64_t nextInRange(
            uint64_t min,
            uint64_t max
        );

        double nextDouble01();
        double nextInRangeDouble(
            double min,
            double max
        );
    };

    //===============================FUNCTIONS=================================

    /*
    * @brief Computes a SplitMix64 hash of the input value
    *
    * @param x Input value to be mixed
    * @return 64 bit hashed value
    */
    static uint64_t splitmix64Hash(
        uint64_t x
    );

    /*
    * @brief Deterministically selects an index for a worker
    *
    * @details
    * Uses a seeded random generator to assign a worker to an entity
    * (bank or firm). The selection depends only on the base
    * seed, the worker's global ID, and the provided stream tag.
    * The assignment must remain identical regardless of system configuration
    *
    * @param baseSeed Base random seed for the simulation.
    * @param workerGlobalId Global ID of the worker being assigned
    * @param countTotal Total number of possible targets to choose from
    * @param streamTag RNG stream identifier
    * @return Selected index between 0 and countTotal - 1
    */
    unsigned int selectForWorker(
        uint64_t baseSeed,
        unsigned int workerGlobalId,
        unsigned int countTotal,
        uint64_t streamTag = 1
    );

    /*
    * @brief Generates a deterministic mixed seed from multiple parameters
    *
    * @details
    * Combines the base seed with run number, timestep, entity ID, and a stream
    * identifier to produce a unique seed. The resulting seed is then passed
    * through a SplitMix64
    *
    * This is to make sure that all random values in the simulation are fully
    * independent of the system configuration
    *
    * @param baseSeed Base random seed of the simulation
    * @param run Current run index
    * @param timestep Current timestep index
    * @param globalId Global ID of the entity associated with the seed
    * @param streamTag RNG stream identifier
    * @return Mixed 64 bit seed value
    */
    static uint64_t makeSeed(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int timestep,
        unsigned int globalId,
        uint64_t streamTag
    );

    /*
    * Generates a percentage value from 0 to 100 based on the mixed seed
    * @note The seed has to be mixed prior to calling this function!
    */
    unsigned short percent0to100FromMixedSeed(
        uint64_t seed
    );

    /*
    * Generates a random value between [minVal, maxVal] based on the mixed seed
    * @note The seed has to be mixed prior to calling this function!
    */
    uint64_t randomInRangeFromSeed(
        uint64_t seed,
        uint64_t minVal,
        uint64_t maxVal
    );

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_utils.cpp                                             *
    **************************************************************************/

    //===============================STRUCTURES================================

    // Debt entry is reused by bank and firm balance sheets
    struct DebtEntry {
        // Global bank ID
        unsigned int lenderBankGlobalId;
        // How much is owed
        uint64_t amount;
    };

    //================================FUNCTIONS================================

    // Returns floor(x * pct / 100) using integer arithmetic
    uint64_t percentFloorU64(
        uint64_t x,
        unsigned short pct
    );

    // Clamps fractions and checks for faulty parameters
    bool errorDetectionAndClamping(
        const unsigned int runs,
        const unsigned int timesteps,
        const uint64_t baseSeed,

        const unsigned int bankCountTotal,
        const unsigned int firmCountTotal,
        const unsigned int workerCountTotal,

        const unsigned int initialBankLiquidity,
        const unsigned int initialFirmLiquidity,
        const unsigned int initialProductionCost,
        const unsigned int wage,

        unsigned short& wageConsumptionPercent,
        unsigned short& profitMultiplierMin,
        unsigned short& profitMultiplierMax,
        unsigned short& shockMultiplierMin,
        unsigned short& shockMultiplierMax,
        unsigned short& minInterestRate,
        unsigned short& maxInterestRate,
        unsigned short& firmShockMinPercent,
        unsigned short& firmShockMaxPercent,
        unsigned short& firmProfitMinPercent,
        unsigned short& firmProfitMaxPercent,
        unsigned short& firmOperatingCostMinPercent,
        unsigned short& firmOperatingCostMaxPercent,
        unsigned short& firmShockProfitMultiplierMinPercent,
        unsigned short& firmShockProfitMultiplierMaxPercent,
        unsigned short& firmInitialLiquidityMultiplierMinPercent,
        unsigned short& firmInitialLiquidityMultiplierMaxPercent,
        unsigned short& firmRareEventProbabilityPercent,
        unsigned short& firmRareEventImpactMinPercent,
        unsigned short& firmRareEventImpactMaxPercent,

        unsigned short& interbankDensity,
        const unsigned short maxInterbankLenderSamplingK,
        unsigned short& maxInterbankLoanPercent,
        unsigned short& maxFirmLoanPercent,
        unsigned short& firmLenderDegree,
        unsigned short& firmRepayPercent,
        unsigned short& bankRepayPercent,
        unsigned int interventionDelay,

        unsigned short& bailInCoveragePercent,
        unsigned short& bailOutCoveragePercent
    );

    // Finds range and starting index for the given rank
    pair<unsigned int, unsigned int> computeRange(
        unsigned int totalCount,
        unsigned int rank,
        unsigned int totalRanks
    );

    // Returns the MPI rank that owns the entity with the given global ID
    unsigned int ownerRankFromGlobalId(
        unsigned int countTotal,
        unsigned int totalRanks,
        unsigned int globalId
    );

    // Generates an interest rate for a bank
    // Used to avoid communication to receive interest rates
    // The interest rate gets regenerated using same parameters as
    // it was originally generated from
    unsigned int getInterestRate(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int globalId,
        uint64_t minVal,
        uint64_t maxVal
    );

    // Makes seed and gets a random percentage value (0 - 100)
    unsigned short randomPercentFromStream(
        const uint64_t baseSeed,
        const unsigned int run,
        const unsigned int step,
        const unsigned int globalId,
        const unsigned int streamTag,
        const unsigned short minPercent,
        const unsigned short maxPercent
    );

    // Returns the total debt amount by adding all entries together
    static uint64_t sumDebtU64(
        const vector<BailInBailOut::DebtEntry>& debts
    );

    // Creates a committed, resized MPI struct type from field descriptors.
    // All block lengths are 1. Frees the raw type internally.
    static MPI_Datatype makeStructType(
        int count,
        const MPI_Aint* offsets,
        const MPI_Datatype* types,
        MPI_Aint extent
    );

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_firms.cpp                                             *
    **************************************************************************/

    //================================STRUCTURES===============================

    struct FirmLoanRequest {
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountRequested;
        unsigned short lenderInterestRate;
        unsigned short pad16;
        unsigned int pad32;
    };

    // Message used to send aggregated workforce cost contributions to the rank
    // that owns the target firm
    struct a0FirmWorkforceDeltaMessage {
        unsigned int firmGlobalId;
        uint64_t delta;
    };

    //================================DATATYPES================================

    MPI_Datatype a0MakeFirmWorkforceDeltaType();

    //================================FUNCTIONS================================

    // Builds the fixed bank neighbor list for a firm
    // The neighbor set is generated and written to outNeighbors
    void buildFirmNeighbors(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int firmGlobalId,
        unsigned int bankCountTotal,
        unsigned int firmLenderDegree,
        uint64_t stream,
        vector<unsigned int>& outNeighbors
    );

    // Computes the random values used in phase A all at once per firm
    // Called at each timestep to generate new values
    void a0BuildFirmPhaseAParameters(
        const uint64_t baseSeed,
        const unsigned int run,
        const unsigned int step,
        const unsigned int firmGlobalStartIndex,
        const unsigned int firmCountForRank,

        const unsigned short firmShockMinPercent,
        const unsigned short firmShockMaxPercent,
        const unsigned short firmProfitMinPercent,
        const unsigned short firmProfitMaxPercent,
        const unsigned short firmOperatingCostMinPercent,
        const unsigned short firmOperatingCostMaxPercent,
        const unsigned short firmShockProfitMultiplierMinPercent,
        const unsigned short firmShockProfitMultiplierMaxPercent,
        const unsigned short firmInitialLiquidityMultiplierMinPercent,
        const unsigned short firmInitialLiquidityMultiplierMaxPercent,
        const unsigned short firmRareEventImpactMinPercent,
        const unsigned short firmRareEventImpactMaxPercent,

        const uint64_t streamFirmShockPercent,
        const uint64_t streamFirmProfitPercent,
        const uint64_t streamFirmOperatingCostPercent,
        const uint64_t streamFirmShockProfitMultiplierPercent,
        const uint64_t streamFirmInitialLiquidityMultiplierPercent,
        const uint64_t streamFirmRareEventProbability,
        const uint64_t streamFirmRareEventImpactPercent,
        const uint64_t streamFirmRareEventSign,

        vector<unsigned short>& localFirmShockMinPercent,
        vector<unsigned short>& localFirmShockMaxPercent,
        vector<unsigned short>& localFirmProfitMinPercent,
        vector<unsigned short>& localFirmProfitMaxPercent,
        vector<unsigned short>& localFirmOperatingCostPercent,
        vector<unsigned short>& localFirmShockProfitMultiplierPercent,
        vector<unsigned short>& localFirmInitialLiquidityMultiplierPercent,
        vector<unsigned short>& localFirmRareEventRollPercent,
        vector<unsigned short>& localFirmRareEventImpactPercent,
        vector<bool>& localFirmRareEventIsPositive
    );

    // Repays firm's outstanding bank debts proportionally based on configured
    // repayment percentage and available liquidity.
    // Records incoming repayments per lender bank and appends each touched
    // bank ID to dirtyBankRepayIds (first write per ID only) using dirty-ID clearing.
    void repayFirmLoansProRata(
        int64_t& firmLiquidity,
        vector<DebtEntry>& debts,
        unsigned short firmRepayPercent,
        vector<int64_t>& bankIncomingFromFirmRepay,
        vector<uint32_t>& dirtyBankRepayIds
    );

    // Computes the total workforce cost owed by each local firm
    // Routes combined wage deltas to the ranks that own those firms
    void a0ComputeFirmWorkforceCost(
        unsigned int firmCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int firmGlobalStartIndex,
        unsigned int firmCountForRank,
        const vector<uint32_t>& localWorkerFirmID,
        unsigned int wage,
        vector<uint64_t>& localFirmWorkforceCost
    );

    // Creates firm loan requests for firms with negative liquidity
    // Routes each request to the rank that owns lender bank
    void a6RequestFirmLoans(
        unsigned int run,
        uint64_t baseSeed,
        unsigned int bankCountTotal,
        unsigned int mpiSize,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        unsigned int firmGlobalStartIndex,
        vector<int64_t>& localFirmLiquidity,
        vector<vector<unsigned int>>& localFirmNeighbors,
        vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank
    );

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_banks.cpp                                             *
    **************************************************************************/

    //===============================STRUCTURES================================

    // Message used to apply a signed liquidity change to a bank
    // Used when routing worker deposits and firm repayments to bank owners.
    struct c1BankDeltaMessage {
        unsigned int bankGlobalId;
        unsigned int pad;
        int64_t delta;
    };

    // Acceptance message sent from lender bank back to firm that requested loan
    // All IDs are global
    struct d2FirmLoanAcceptance {
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };

    // Message representing an interbank loan repayment sent to the lender bank
    // The amount is credited to the bank
    struct e1InterbankRepaymentMessage {
        unsigned int lenderBankGlobalId;
        unsigned int pad;
        uint64_t amount;
    };

    // Request sent by bank to potential lender bank for interbank liquidity
    // Includes the requested amount and the lender's interest rate
    struct e3InterbankLoanRequest {
        unsigned int borrowerBankGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountRequested;
        unsigned short lenderInterestRate;
        unsigned short pad16;
        unsigned int pad32;
    };

    // Acceptance message sent by lender bank after granting interbank loan
    // Borrower applies granted amount and records the new debt locally
    struct e3InterbankLoanAcceptance {
        unsigned int borrowerBankGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };




    // Local helper record representing one possible interbank lender
    // Stores lender ID, interest rate, and current liquidity for ranking
    struct e3LenderOption {
        unsigned int lenderGlobalId;
        unsigned short interestRate;
        int64_t liquidity;
    };

    //================================DATATYPES================================

    MPI_Datatype d1MakeFirmLoanRequestType();
    MPI_Datatype d4MakeFirmLoanAcceptanceType();
    MPI_Datatype c1MakeBankDeltaMessageType();
    MPI_Datatype e1MakeInterbankRepaymentType();
    MPI_Datatype e3MakeInterbankLoanRequestType();
    MPI_Datatype e3MakeInterbankLoanAcceptanceType();

    //================================FUNCTIONS================================

    // Transfers worker deposits and firm loan repayments to owning bank ranks
    // and applies the resulting liquidity changes to local banks.
    // wageDepositBuffer and wageNonZeroIds are pre-computed per run and never
    // modified here. bankIncomingFromFirmRepay is zeroed only for dirty IDs via dirtyBankRepayIds.
    // c1ToRank is a pre-allocated per-policy bucket of size mpiSize.
    void c1SendMoneyToBanks(
        unsigned int bankCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        vector<int64_t>& localBankLiquidity,
        const vector<uint32_t>& wageNonZeroIds,
        const vector<int64_t>& wageDepositBuffer,
        vector<int64_t>& bankIncomingFromFirmRepay,
        vector<uint32_t>& dirtyBankRepayIds,
        vector<vector<c1BankDeltaMessage>>& c1ToRank
    );

    // Sends firm loan requests to the ranks that own the target banks and
    // collects all incoming requests that are processed locally
    void d1SendFirmLoanRequests(
        unsigned int mpiRank,
        unsigned int mpiSize,
        vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank,
        vector<FirmLoanRequest>& receivedLoanRequests
    );

    // Processes incoming firm loan requests for banks owned by this rank,
    // grants loans as possible (liquidity limits)
    void d2ProcessFirmLoanRequests(
        unsigned int firmCountTotal,
        unsigned int bankCountForRank,
        unsigned int bankGlobalStartIndex,
        unsigned int mpiSize,
        unsigned short maxFirmLoanPercent,
        vector<int64_t>& localBankLiquidity,
        vector<FirmLoanRequest>& receivedLoanRequests,
        vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
        vector<uint64_t>& bankFirmLoanOutflow
    );

    // Sends loan acceptance messages from lender banks back to the ranks
    // that own the borrowing firms
    void d4SendFirmLoanAcceptances(
        unsigned int mpiRank,
        unsigned int mpiSize,
        vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
        vector<d2FirmLoanAcceptance>& receivedAcceptances
    );

    // Applies accepted firm loans by increasing firm liquidity and recording
    // the new debt entries
    void d5ApplyFirmLoanAcceptances(
        unsigned int firmGlobalStartIndex,
        vector<int64_t>& localFirmLiquidity,
        vector<vector<DebtEntry>>& localFirmDebts,
        vector<d2FirmLoanAcceptance>& receivedAcceptances
    );

    // Repays outstandi ng interbank debts proportionally using available bank
    // liquidity and routes repayment messages to lender banks
    // e1ToRank is a pre-allocated per-policy bucket of size mpiSize.
    void e1RepayInterbankLoans(
        unsigned int bankCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        unsigned short bankRepayPercent,
        vector<int64_t>& localBankLiquidity,
        vector<vector<DebtEntry>>& localBankDebts,
        vector<vector<e1InterbankRepaymentMessage>>& e1ToRank
    );

    // Applies interest to all existing debts using the regenerated interest
    // rates per bank
    void e2ApplyInterestOnAllLoans(
        vector<vector<DebtEntry>>& eLocalBankDebts,
        vector<vector<DebtEntry>>& localFirmDebts,
        uint64_t baseSeed,
        unsigned int run,
        uint64_t minVal,
        uint64_t maxVal
    );

    // Banks with negative liquidity borrow from neighboring banks in the
    // interbank graph. E.3.1 uses MPI_Allgatherv to populate lenderLiquidityView
    // with every bank's current liquidity in one collective (replaces the prior
    // neighbor-liquidity alltoallv that scaled as O(mpiSize^2) in message
    // startup cost). bankAllgatherRecvCounts and bankAllgatherDispls are
    // pre-computed per run from computeRange and describe how ranks partition
    // banks (supports uneven splits).
    //
    // K lenders are selected randomly (partial Fisher-Yates) using a seed
    // derived from baseSeed, run, step, and borrower ID.
    // samplingIndices is a pre-allocated identity-init buffer of size
    // `degree`; the function performs at most K swaps per bank and restores
    // the buffer to identity by reverse-swap before returning per bank.
    // reqsByLocalLender is a pre-allocated bucket of size bankCountForRank;
    // each inner vector is cleared at the start of E.3.3 and refilled with
    // received requests for that step.
    // e3RequestsToRank and e3AcceptancesToRank are pre-allocated per-policy
    // buckets of size mpiSize for E.3.3 and E.3.5 outbound messages.
    void e3BorrowInterbankIfNegativeLiquidity(
        unsigned int bankCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        unsigned int step,
        unsigned short maxInterbankLenderSamplingK,
        unsigned short maxInterbankLoanPercent,
        uint64_t baseSeed,
        unsigned int run,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        vector<vector<unsigned int>>& localBankNeighbors,
        const vector<int>& bankAllgatherRecvCounts,
        const vector<int>& bankAllgatherDispls,
        vector<int64_t>& localBankLiquidity,
        vector<vector<DebtEntry>>& localBankDebts,
        vector<int64_t>& lenderLiquidityView,
        vector<unsigned int>& samplingIndices,
        vector<vector<e3InterbankLoanRequest>>& reqsByLocalLender,
        vector<vector<e3InterbankLoanRequest>>& e3RequestsToRank,
        vector<vector<e3InterbankLoanAcceptance>>& e3AcceptancesToRank
    );

    void e3BorrowInterbankIfNegativeLiquidityDense(
        unsigned int bankCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        unsigned int step,
        unsigned short maxInterbankLenderSamplingK,
        unsigned short maxInterbankLoanPercent,
        uint64_t baseSeed,
        unsigned int run,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        vector<vector<unsigned int>>& localBankNeighbors,
        const vector<int>& bankAllgatherRecvCounts,
        const vector<int>& bankAllgatherDispls,
        vector<int64_t>& localBankLiquidity,
        vector<vector<DebtEntry>>& localBankDebts,
        vector<int64_t>& lenderLiquidityView,
        vector<unsigned int>& samplingIndices,
        vector<vector<e3InterbankLoanRequest>>& reqsByLocalLender,
        vector<vector<e3InterbankLoanRequest>>& e3RequestsToRank,
        vector<vector<e3InterbankLoanAcceptance>>& e3AcceptancesToRank
    );


    void f1f2UpdateBankDistressAndApplyIntervention(
        unsigned int bankCountForRank,
        unsigned int interventionDelay,
        unsigned int policy,
        unsigned short bailInCoveragePercent,
        unsigned short bailOutCoveragePercent,
        vector<int64_t>& localBankLiquidity,
        vector<vector<DebtEntry>>& localBankDebts,
        vector<unsigned int>& bankDistressCount,
        uint64_t& bailMoneyThisStep,
        uint64_t& graceMoneyThisStep
    );

    /**************************************************************************
    * This header section defines and implements                              *
    * all templated functions used by the program                             *
    **************************************************************************/

    // Variable-length all-to-all exchange using non-blocking point-to-point.
    // mpiRank must be passed by caller (no internal MPI_Comm_rank call).
    // tag defaults to 0; pass a distinct value only when two concurrent
    // exchanges could otherwise mix messages (rare with sequential phases).
    //
    // The sendCounts / recvCounts / recvOffsets scratch vectors are hoisted
    // to function-local static storage: each template instantiation owns one
    // set, resized lazily on first use. This eliminates 3 * mpiSize int
    // allocations per call, which matters at high rank counts (7 calls/step).
    // MPI is single-threaded per rank in this program, so static is safe.
    template<typename T>
    void allToAllvExchange(
        unsigned int mpiRank,
        unsigned int mpiSize,
        vector<vector<T>>& toRank,
        vector<T>& received,
        MPI_Datatype mpiType,
        int tag = 0
    ) {
        (void)mpiType;
        (void)tag;

        if (mpiSize == 1) {
            if (toRank[0].size() > (
                size_t
            )numeric_limits<int>::max()) {
                cerr << "allToAllvExchange send count exceeds MPI int limit on rank "
                    << mpiRank << " for destination 0\n";
                MPI_Abort(
                    MPI_COMM_WORLD,
                    1
                );
            }

            received.clear();
            received.insert(
                received.end(),
                toRank[0].begin(),
                toRank[0].end()
            );
            toRank[0].clear();
            return;
        }

        static vector<int> sendCounts;
        static vector<int> recvCounts;
        static vector<int> recvOffsets;
        if (sendCounts.size() != mpiSize) sendCounts.resize(
            mpiSize
        );
        if (recvCounts.size() != mpiSize) recvCounts.resize(
            mpiSize
        );
        if (recvOffsets.size() != mpiSize) {
            recvOffsets.resize(
                mpiSize
            );
        }

        for (unsigned int r = 0; r < mpiSize; r++) {
            if (toRank[r].size() > (
                size_t
            )numeric_limits<int>::max()) {
                cerr << "allToAllvExchange send count exceeds MPI int limit on rank "
                    << mpiRank << " for destination " << r << "\n";
                MPI_Abort(
                    MPI_COMM_WORLD,
                    1
                );
            }

            sendCounts[r] = (int)toRank[r].size();
        }

        MPI_Alltoall(
            sendCounts.data(), 1, MPI_INT,
            recvCounts.data(), 1, MPI_INT,
            MPI_COMM_WORLD
        );

        size_t recvTotalSize = 0;
        for (unsigned int r = 0; r < mpiSize; r++) {
            recvOffsets[r] = (int)recvTotalSize;
            recvTotalSize += (size_t)recvCounts[r];

            if (recvTotalSize > (
                size_t
            )numeric_limits<int>::max()) {
                cerr << "allToAllvExchange receive total exceeds MPI int limit on rank "
                    << mpiRank << "\n";
                MPI_Abort(
                    MPI_COMM_WORLD,
                    1
                );
            }
        }

        int recvTotal = (int)recvTotalSize;

        received.clear();
        received.resize(
            (unsigned int)recvTotal
        );

        static vector<MPI_Request> requests;
        requests.clear();
        requests.reserve(
            (size_t)2 * mpiSize
        );

        for (unsigned int r = 0; r < mpiSize; r++) {
            if (recvCounts[r] == 0 || r == mpiRank) {
                continue;
            }
            MPI_Request req;
            MPI_Irecv(
                received.data() + recvOffsets[r],
                recvCounts[r],
                mpiType,
                (int)r,
                tag,
                MPI_COMM_WORLD,
                &req
            );
            requests.push_back(
                req
            );
        }

        for (unsigned int r = 0; r < mpiSize; r++) {
            if (sendCounts[r] == 0 || r == mpiRank) {
                continue;
            }
            MPI_Request req;
            MPI_Isend(
                toRank[r].data(),
                sendCounts[r],
                mpiType,
                (int)r,
                tag,
                MPI_COMM_WORLD,
                &req
            );
            requests.push_back(
                req
            );
        }

        if (sendCounts[mpiRank] > 0) {
            int selfOffset = recvOffsets[mpiRank];
            for (int i = 0; i < sendCounts[mpiRank]; i++) {
                received[(size_t)selfOffset + (size_t)i] =
                    toRank[mpiRank][(size_t)i];
            }
        }

        if (!requests.empty()) {
            MPI_Waitall(
                (int)requests.size(),
                requests.data(),
                MPI_STATUSES_IGNORE
            );
        }

        for (unsigned int r = 0; r < mpiSize; r++) {
            toRank[r].clear();
        }
    }


};

#endif
