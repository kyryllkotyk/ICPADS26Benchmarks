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

#include "pm2_bailinbailout.h"

// Run bail-in and bail-out policies back to back with same seed inputs
void BailInBailOut::runBenchmarkInAndOut(
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
) {

    // Validate and clamp input parameters before allocating state

    if (!errorDetectionAndClamping(
        runs,
        timesteps,
        baseSeed,
        bankCountTotal,
        firmCountTotal,
        workerCountTotal,
        initialBankLiquidity,
        initialFirmLiquidity,
        initialProductionCost,
        wage,
        wageConsumptionPercent,
        profitMultiplierMin,
        profitMultiplierMax,
        shockMultiplierMin,
        shockMultiplierMax,
        minInterestRate,
        maxInterestRate,
        firmShockMinPercent,
        firmShockMaxPercent,
        firmProfitMinPercent,
        firmProfitMaxPercent,
        firmOperatingCostMinPercent,
        firmOperatingCostMaxPercent,
        firmShockProfitMultiplierMinPercent,
        firmShockProfitMultiplierMaxPercent,
        firmInitialLiquidityMultiplierMinPercent,
        firmInitialLiquidityMultiplierMaxPercent,
        firmRareEventProbabilityPercent,
        firmRareEventImpactMinPercent,
        firmRareEventImpactMaxPercent,
        interbankDensity,
        maxInterbankLenderSamplingK,
        maxInterbankLoanPercent,
        maxFirmLoanPercent,
        firmLenderDegree,
        firmRepayPercent,
        bankRepayPercent,
        interventionDelay,
        bailInCoveragePercent,
        bailOutCoveragePercent
    )) {
        return;
    }

    // Get MPI rank information used for ownership and reductions

    int mpiRank;
    MPI_Comm_rank(
        MPI_COMM_WORLD,
        &mpiRank
    );

    int mpiSize;
    MPI_Comm_size(
        MPI_COMM_WORLD,
        &mpiSize
    );


    // Split bank, firm, and worker ranges across ranks

    // Banks
    unsigned int bankCountForRank, bankGlobalStartIndex;
    // Each rank owns contiguous ranges for all 3 entity types
    pair<unsigned int, unsigned int> bankPair = computeRange(
        bankCountTotal,
        mpiRank,
        mpiSize
    );
    bankCountForRank = bankPair.first;
    bankGlobalStartIndex = bankPair.second;

    // Firms
    unsigned int firmCountForRank, firmGlobalStartIndex;
    pair<unsigned int, unsigned int> firmPair = computeRange(
        firmCountTotal,
        mpiRank,
        mpiSize
    );
    firmCountForRank = firmPair.first;
    firmGlobalStartIndex = firmPair.second;

    // Workers
    unsigned int workerCountForRank, workerGlobalStartIndex;
    pair<unsigned int, unsigned int> workerPair = computeRange(
        workerCountTotal,
        mpiRank,
        mpiSize
    );
    workerCountForRank = workerPair.first;
    workerGlobalStartIndex = workerPair.second;
    for (int run = 0; run < runs; run++) {
        uint64_t totalBailOutMoney = 0, totalBailInMoney = 0;

        // Per-run buffers reused by both policy modes
        // These buffers depend only on (baseSeed, run) and worker/bank
        // assignments; they do NOT depend on policy. Built once per run
        // and reused across both Bail-In and Bail-Out iterations.

        vector<uint32_t> localWorkerBankID(
            workerCountForRank
        );
        vector<uint32_t> localWorkerFirmID(
            workerCountForRank
        );
        for (int i = 0; i < (int)localWorkerBankID.size(); i++) {
            localWorkerBankID[i] = selectForWorker(
                baseSeed,
                i + workerGlobalStartIndex,
                bankCountTotal,
                STREAM_WORKER_BANK_ASSIGN
            );

            localWorkerFirmID[i] = selectForWorker(
                baseSeed,
                i + workerGlobalStartIndex,
                firmCountTotal,
                STREAM_WORKER_FIRM_ASSIGN
            );
        }

        unsigned int degree = (unsigned int)(((uint64_t)(bankCountTotal
            - 1) * interbankDensity) / 100);

        vector<vector<unsigned int>> localBankNeighbors(
            bankCountForRank
        );
        vector<uint32_t> usedBankStamp(
            bankCountTotal,
            0
        );
        uint32_t currentStamp = 1;

        for (unsigned int localBorrower = 0; localBorrower <
            bankCountForRank; localBorrower++) {

            unsigned int borrowerGlobalId = bankGlobalStartIndex
                + localBorrower;

            vector<unsigned int>& neighbors = localBankNeighbors[localBorrower];
            neighbors.clear();

            if (degree == 0) {
                continue;
            }

            neighbors.reserve(
                degree
            );

            uint64_t state = makeSeed(
                baseSeed,
                run,
                0,
                borrowerGlobalId,
                STREAM_BANK_GRAPH_BUILD
            );

            if (currentStamp == 0) {
                fill(
                    usedBankStamp.begin(),
                    usedBankStamp.end(),
                    0
                );
                currentStamp = 1;
            }

            uint32_t stamp = currentStamp++;
            usedBankStamp[borrowerGlobalId] = stamp;

            while (neighbors.size() < degree) {
                state = splitmix64Hash(
                    state
                );
                unsigned int candidate = (unsigned int)(state %
                    bankCountTotal);

                if (usedBankStamp[candidate] != stamp) {
                    usedBankStamp[candidate] = stamp;
                    neighbors.push_back(
                        candidate
                    );
                }
            }
        }

        vector<vector<unsigned int>> localFirmNeighbors(
            firmCountForRank
        );
        for (unsigned int localFirm = 0; localFirm < firmCountForRank; localFirm++) {
            unsigned int firmGlobalId = firmGlobalStartIndex + localFirm;
            buildFirmNeighbors(
                baseSeed,
                (unsigned int)run,
                firmGlobalId,
                bankCountTotal,
                (unsigned int)firmLenderDegree,
                STREAM_FIRM_BANK_GRAPH_BUILD,
                localFirmNeighbors[localFirm]
            );
        }

        vector<uint64_t> localFirmWorkforceCost(
            firmCountForRank,
            0
        );
        a0ComputeFirmWorkforceCost(
            firmCountTotal,
            mpiRank,
            mpiSize,
            firmGlobalStartIndex,
            firmCountForRank,
            localWorkerFirmID,
            wage,
            localFirmWorkforceCost
        );

        vector<int64_t> wageDepositBuffer(
            bankCountTotal,
            0
        );
        vector<uint32_t> wageNonZeroIds;
        wageNonZeroIds.reserve(
            bankCountTotal
        );
        {
            uint64_t depositPerW = ((uint64_t)wage *
                (uint64_t)(100 - wageConsumptionPercent)) / 100ULL;
            for (unsigned int wID = 0; wID < workerCountForRank; wID++) {
                uint32_t bankGID = (uint32_t)localWorkerBankID[wID];
                if (wageDepositBuffer[bankGID] == 0) wageNonZeroIds.push_back(
                    bankGID
                );
                wageDepositBuffer[bankGID] += (int64_t)depositPerW;
            }
        }

        // E.3.1 now uses MPI_Allgatherv to gather every bank's liquidity onto
        // every rank. Pre-compute the recv counts and displacements once per
        // run (depends only on bankCountTotal and mpiSize via computeRange).
        // These describe how banks are partitioned: ranks [0, remainder) own
        // (base+1) banks each, the rest own `base` banks. Displacements are
        // the global start index of each rank's bank range.
        // Bank allgather layout mirrors computeRange ownership
    vector<int> bankAllgatherRecvCounts(
            (unsigned int)mpiSize,
            0
        );
        vector<int> bankAllgatherDispls(
            (unsigned int)mpiSize,
            0
        );
        {
            int runningDispl = 0;
            for (unsigned int r = 0; r < (unsigned int)mpiSize; r++) {
                pair<unsigned int, unsigned int> rr = computeRange(
                    bankCountTotal,
                    r,
                    (unsigned int)mpiSize
                );
                bankAllgatherRecvCounts[r] = (int)rr.first;
                bankAllgatherDispls[r] = runningDispl;
                runningDispl += (int)rr.first;
            }
        }

        // Policy = 0 -> Bail-In
        // Policy = 1 -> Bail-Out
        for (int policy = 0; policy <= 1; policy++) {
            uint64_t totalDeficit = 0;
            // Per-policy local entity state
            // Bank Arrays
            // Local bank and firm liquidity start fresh for each policy
            vector<int64_t> localBankLiquidity(
                bankCountForRank,
                initialBankLiquidity
            );
            vector<unsigned short> localBankInterestRate(
                bankCountForRank
            );
            vector<FirmLoanRequest> receivedLoanRequests(
                0
            );

            vector<vector<d2FirmLoanAcceptance>> firmLoanAcceptancesToRank(
                mpiSize
            );
            vector<uint64_t> bankFirmLoanOutflow(
                bankCountForRank
            );
            vector<d2FirmLoanAcceptance> receivedAcceptances;
            vector<vector<DebtEntry>> eLocalBankDebts(
                bankCountForRank
            );
            vector<unsigned int> bankDistressCount(
                bankCountForRank,
                0
            );
            vector<int64_t> bankIncomingFromFirmRepay(
                bankCountTotal,
                0
            );

            // Firm Arrays
            vector<unsigned short> localFirmShockMinPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmShockMaxPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmProfitMinPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmProfitMaxPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmOperatingCostPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmShockProfitMultiplierPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmInitialLiquidityMultiplierPercent(
                firmCountForRank, 0);
            vector<unsigned short> localFirmRareEventImpactPercent(
                firmCountForRank, 0);
            vector<int64_t> localFirmLiquidity(
                firmCountForRank,
                initialFirmLiquidity
            );
            vector<uint64_t> localFirmProductionCost(
                firmCountForRank,
                initialProductionCost
            );

            for (unsigned int localFirm = 0; localFirm < firmCountForRank; localFirm++) {
                unsigned int firmGlobalId = firmGlobalStartIndex + localFirm;
                unsigned short initialLiquidityMultiplier = randomPercentFromStream(
                    baseSeed,
                    (unsigned int)run,
                    0,
                    firmGlobalId,
                    STREAM_FIRM_INITIAL_LIQUIDITY_MULTIPLIER_PERCENT,
                    firmInitialLiquidityMultiplierMinPercent,
                    firmInitialLiquidityMultiplierMaxPercent
                );

                localFirmLiquidity[localFirm] =
                    (int64_t)(((uint64_t)initialFirmLiquidity *
                        (uint64_t)initialLiquidityMultiplier) / 100ULL);
            }

            vector<unsigned short> localFirmRareEventRollPercent(
                firmCountForRank
            );
            vector<bool> localFirmRareEventIsPositive(
                firmCountForRank
            );
            vector<vector<DebtEntry>> localFirmDebts(
                firmCountForRank
            );
            vector<vector<FirmLoanRequest>> firmLoanRequestsToRank(
                mpiSize
            );

            // (Worker arrays, firm/bank neighbor graphs, workforce cost,
            //  wageDepositBuffer/wageNonZeroIds, and the bank allgather
            //  routing vectors (bankAllgatherRecvCounts/bankAllgatherDispls)
            //  are hoisted to per-run scope above the policy loop.)

            // Initialize bank liquidity, rates, and interbank graph state

            for (unsigned int i = 0; i < bankCountForRank; i++) {
                // Get the global bank ID
                unsigned int bankGlobalId = bankGlobalStartIndex + i;

                // Find the seed for this bank's interest rate
                uint64_t bankSeed = makeSeed(
                    baseSeed,
                    run,
                    0,
                    bankGlobalId,
                    STREAM_INTEREST_RATE_SELECTION
                );

                localBankInterestRate[i] = randomInRangeFromSeed(
                    bankSeed,
                    minInterestRate,
                    maxInterestRate
                );
            }

            // Per-policy scratch buffers reused across timesteps
            // localBankNeighbors / localFirmNeighbors / localFirmWorkforceCost
            // / wageDepositBuffer / wageNonZeroIds / bankAllgatherRecvCounts
            // / bankAllgatherDispls / degree are all hoisted to per-run scope
            // above. The buffers below hold per-policy mutable scratch state
            // and live for the timestep loop.

            // Size bankCountTotal so MPI_Allgatherv fills it in one shot.
            // Allgather populates every slot, so no sentinel init is needed.
            vector<int64_t> lenderLiquidityView(
                bankCountTotal,
                0
            );
            vector<uint32_t> dirtyBankRepayIds;
            dirtyBankRepayIds.reserve(
                bankCountForRank
            );

            // E.3 K-sampling: identity-init once; per-bank Fisher-Yates
            // restores via reverse swaps so the buffer is reusable.
            vector<unsigned int> e3SamplingIndices(
                degree
            );
            for (unsigned int i = 0; i < degree; i++) {
                e3SamplingIndices[i] = i;
            }

            // E.3.3 lender-side request bucket: persistent buffer cleared
            // each step. Inner vectors retain capacity across timesteps.
            vector<vector<e3InterbankLoanRequest>> reqsByLocalLender(
                bankCountForRank
            );

            // Hoisted outbound rank buckets. allToAllvExchange clears the
            // inner vectors at the end of each call, preserving capacity; the
            // outer allocation of size mpiSize happens once per policy. These
            // replace the prior per-step `vector<vector<T>>(mpiSize)` inside
            // c1 / e1 / e3 which allocated mpiSize inner vectors every step.
            vector<vector<c1BankDeltaMessage>> c1ToRank(
                mpiSize
            );
            vector<vector<e1InterbankRepaymentMessage>> e1ToRank(
                mpiSize
            );
            vector<vector<e3InterbankLoanRequest>> e3RequestsToRank(
                mpiSize
            );
            vector<vector<e3InterbankLoanAcceptance>> e3AcceptancesToRank(
                mpiSize
            );

            // Main timestep loop for all simulation phases

            auto timestepStart = chrono::high_resolution_clock::now();

            for (int step = 0; step < timesteps; step++) {

                // Phase A updates firm production and worker wages

                // At the start of the phase, calculate all random values
                a0BuildFirmPhaseAParameters(
                    baseSeed,
                    (unsigned int)run,
                    (unsigned int)step,
                    firmGlobalStartIndex,
                    firmCountForRank,

                    firmShockMinPercent,
                    firmShockMaxPercent,
                    firmProfitMinPercent,
                    firmProfitMaxPercent,
                    firmOperatingCostMinPercent,
                    firmOperatingCostMaxPercent,
                    firmShockProfitMultiplierMinPercent,
                    firmShockProfitMultiplierMaxPercent,
                    firmInitialLiquidityMultiplierMinPercent,
                    firmInitialLiquidityMultiplierMaxPercent,
                    firmRareEventImpactMinPercent,
                    firmRareEventImpactMaxPercent,

                    STREAM_FIRM_SHOCK_PERCENT,
                    STREAM_FIRM_PROFIT_PERCENT,
                    STREAM_FIRM_OPERATING_COST_PERCENT,
                    STREAM_FIRM_SHOCK_PROFIT_MULTIPLIER_PERCENT,
                    STREAM_FIRM_INITIAL_LIQUIDITY_MULTIPLIER_PERCENT,
                    STREAM_FIRM_RARE_EVENT_PROBABILITY,
                    STREAM_FIRM_RARE_EVENT_IMPACT_PERCENT,
                    STREAM_FIRM_RARE_EVENT_SIGN,

                    localFirmShockMinPercent,
                    localFirmShockMaxPercent,
                    localFirmProfitMinPercent,
                    localFirmProfitMaxPercent,
                    localFirmOperatingCostPercent,
                    localFirmShockProfitMultiplierPercent,
                    localFirmInitialLiquidityMultiplierPercent,
                    localFirmRareEventRollPercent,
                    localFirmRareEventImpactPercent,
                    localFirmRareEventIsPositive
                );

                for (int f = 0; f < firmCountForRank; f++) {
                    unsigned int firmGlobalId = f + firmGlobalStartIndex;

                    // A.1: Find the random shock value
                    uint64_t shockSignSeed = makeSeed(
                        baseSeed,
                        run,
                        step,
                        firmGlobalId,
                        STREAM_SHOCK_SIGN
                    );

                    short shock = (short)randomPercentFromStream(
                        baseSeed,
                        run,
                        step,
                        firmGlobalId,
                        STREAM_SHOCK_GENERATION,
                        localFirmShockMinPercent[f],
                        localFirmShockMaxPercent[f]
                    );

                    if (shockSignSeed & 1ULL) {
                        shock *= -1;
                    }

                    shock = (short)(((int)shock *
                        (int)localFirmShockProfitMultiplierPercent[f]) / 100);

                    // A.2: Affect production cost using shock value.
                    // Clamp below zero so negative shocks cannot wrap into a huge uint64_t.
                    int productionCostPercent = 100 + (int)shock;
                    if (productionCostPercent < 0) {
                        productionCostPercent = 0;
                    }

                    localFirmProductionCost[f] =
                        (localFirmProductionCost[f] * (uint64_t)productionCostPercent) / 100ULL;

                    // A.3: Find the random profit multiplier and affect profit
                    short profitMultiplier = (short)randomPercentFromStream(
                        baseSeed,
                        run,
                        step,
                        firmGlobalId,
                        STREAM_PROFIT_MULTIPLIER_GENERATION,
                        localFirmProfitMinPercent[f],
                        localFirmProfitMaxPercent[f]
                    );

                    profitMultiplier = (short)(((int)profitMultiplier *
                        (int)localFirmShockProfitMultiplierPercent[f]) / 100);

                    int64_t operatingCost = (int64_t)percentFloorU64(
                        localFirmProductionCost[f],
                        localFirmOperatingCostPercent[f]
                    );

                    // A.4: Apply profit and subtract workforce + operating costs
                    uint64_t productionRevenue =
                        (localFirmProductionCost[f] * (uint64_t)profitMultiplier) / 100ULL;

                    localFirmLiquidity[f] +=
                        (int64_t)productionRevenue -
                        (int64_t)localFirmWorkforceCost[f] -
                        operatingCost;

                    // A.4.5: Apply rare event if triggered this timestep
                    if (localFirmRareEventRollPercent[f] < firmRareEventProbabilityPercent) {
                        int64_t rareEventImpact = (int64_t)percentFloorU64(
                            localFirmProductionCost[f],
                            localFirmRareEventImpactPercent[f]
                        );

                        if (localFirmRareEventIsPositive[f]) {
                            localFirmLiquidity[f] += rareEventImpact;
                        }
                        else {
                            localFirmLiquidity[f] -= rareEventImpact;
                        }
                    }

                    // A.5: Repay bank loans (firmRepayPercent % of the loan)
                    // If liquidity is smaller than the percentage of the loan
                    // -> Pay all of liquidity
                    // If liquidity is negative
                    // -> Don't pay at all
                    repayFirmLoansProRata(
                        localFirmLiquidity[f],
                        localFirmDebts[f],
                        firmRepayPercent,
                        bankIncomingFromFirmRepay,
                        dirtyBankRepayIds
                    );
                }

                // A.5.5: Bank pays employees their wage.
                // They keep a percentage of it in their bank accounts
                // A bank employee deposits into the bank they work for
                for (int bankID = 0; bankID < bankCountForRank; bankID++) {
                    localBankLiquidity[bankID] -=
                        (int64_t)bankEmployeeWage *
                        (int64_t)bankWorkerCount *
                        (int64_t)wageConsumptionPercent / 100;
                }

                // A.6: Request loans if needed

                // Clear request buffer
                for (int r = 0; r < mpiSize; r++) {
                    firmLoanRequestsToRank[r].clear();
                }

                // Firms with negative liquidity request loans from selected banks
                a6RequestFirmLoans(
                    run,
                    baseSeed,
                    bankCountTotal,
                    mpiSize,
                    minInterestRate,
                    maxInterestRate,
                    firmGlobalStartIndex,
                    localFirmLiquidity,
                    localFirmNeighbors,
                    firmLoanRequestsToRank
                );

                // Phase C sends deposits and repayments to bank owners

                // Route worker deposits and firm repayments to bank owners
                c1SendMoneyToBanks(
                    bankCountTotal,
                    (unsigned int)mpiRank,
                    (unsigned int)mpiSize,
                    bankGlobalStartIndex,
                    bankCountForRank,
                    localBankLiquidity,
                    wageNonZeroIds,
                    wageDepositBuffer,
                    bankIncomingFromFirmRepay,
                    dirtyBankRepayIds,
                    c1ToRank
                );

                // Phase D processes firm loan requests and acceptances

                // Clear before starting the requesting process
                receivedLoanRequests.clear();

                // D.1: Send loan request buffers to the banks
                // Send firm loan requests to lender bank owner ranks
                d1SendFirmLoanRequests(
                    mpiRank,
                    mpiSize,
                    firmLoanRequestsToRank,
                    receivedLoanRequests
                );

                // D.2: Banks process loan requests by ID
                // Lender banks grant requests up to available lending capacity
                d2ProcessFirmLoanRequests(
                    firmCountTotal,
                    bankCountForRank,
                    bankGlobalStartIndex,
                    mpiSize,
                    maxFirmLoanPercent,
                    localBankLiquidity,
                    receivedLoanRequests,
                    firmLoanAcceptancesToRank,
                    bankFirmLoanOutflow
                );

                // D.3: Banks update their liquidity to reflect accepted loans
                for (unsigned int b = 0; b < bankCountForRank; b++) {
                    localBankLiquidity[b] -= (int64_t)bankFirmLoanOutflow[b];

                    // Liquidity should not go negative from lending
                    if (bankFirmLoanOutflow[b] > 0 && localBankLiquidity[b] < 0) {
                        printf(
                            "ERROR: Bank liquidity negative after D3\n"
                        );
                    }
                }

                // D.4: Banks send back acceptances
                // Return loan acceptance messages to borrower firm owners
                d4SendFirmLoanAcceptances(
                    mpiRank,
                    mpiSize,
                    firmLoanAcceptancesToRank,
                    receivedAcceptances
                );

                // D.5: If accepted, firms update liquidity and loan list
                // Borrower firms apply granted cash and record debt entries
                d5ApplyFirmLoanAcceptances(
                    firmGlobalStartIndex,
                    localFirmLiquidity,
                    localFirmDebts,
                    receivedAcceptances
                );

                // Phase E updates interbank debts and borrowing

                // E.1: Repay interbank loans (bankRepayPercent % of them)
                // The repayments are saved in a buffer
                // The subtraction is applied instantly

                // Banks repay outstanding interbank debts before new borrowing
                e1RepayInterbankLoans(
                    bankCountTotal,
                    mpiRank,
                    mpiSize,
                    bankGlobalStartIndex,
                    bankCountForRank,
                    bankRepayPercent,
                    localBankLiquidity,
                    eLocalBankDebts,
                    e1ToRank
                );

                // E.2: Accrue interest on outstanding loans
                // Applied by the borrower
                // Iterate debt lists for firms and banks, and update loan amount
                // Apply deterministic interest growth to firm and bank debts
                e2ApplyInterestOnAllLoans(
                    eLocalBankDebts,
                    localFirmDebts,
                    baseSeed,
                    run,
                    minInterestRate,
                    maxInterestRate
                );

                // E.3: If the bank's liquidity is negative, attempt to
                // borrow from another bank (by sending a request message)
                // Negative-liquidity banks borrow from interbank neighbors
                e3BorrowInterbankIfNegativeLiquidity(
                    bankCountTotal,
                    mpiRank,
                    mpiSize,
                    bankGlobalStartIndex,
                    bankCountForRank,
                    (unsigned int)step,
                    maxInterbankLenderSamplingK,
                    maxInterbankLoanPercent,
                    baseSeed,
                    (unsigned int)run,
                    minInterestRate,
                    maxInterestRate,
                    localBankNeighbors,
                    bankAllgatherRecvCounts,
                    bankAllgatherDispls,
                    localBankLiquidity,
                    eLocalBankDebts,
                    lenderLiquidityView,
                    e3SamplingIndices,
                    reqsByLocalLender,
                    e3RequestsToRank,
                    e3AcceptancesToRank
                );

                // Phase F applies distress detection and intervention policy

                // F.1: Insolvency check for each local bank
                // Insolvent means bank's liquidity is smaller than the total debt
                // If insolvent, increment the bank's distress counter
                // If solvent, reset the distress counter to 0
                // F.2: If insolvent and the bank has been distressed for at least
                // interventionDelay timesteps, apply policy (once per distress episode)
                //
                // Policy 0 (Bail-In):
                //  Deficit = amount needed
                //  L = total interbank liabilities (debts)
                //  haircutRatio = min(1.0, deficit / L)
                //  haircutFraction = haircutRatio * bailInCoveragePercent / 100.0
                //
                //
                // Policy 1 (Bail-Out, interbank only):
                //  Each interbank loan owed by this bank is reduced multiplicatively
                //  by (1 - haircutFraction).
                //  Deficit = amount needed to restore solvency
                //  Injection = Deficit * bailOutCoveragePercent / 100.0
                //
                //  The injection is added directly to the bank's liquidity.
                //
                // After the policy is applied, mark the intervention as applied for
                // this distress episode so it is not applied again until the bank
                // returns to solvency.
                //
                // For both policies, track how much bail was granted this step
                // (cash injected for bail-out or total debt removed for bail-in)
                uint64_t bailMoneyThisStep = 0, deficitThisStep = 0;
                // Update distress counters and apply selected policy if due
                f1f2UpdateBankDistressAndApplyIntervention(
                    bankCountForRank,
                    interventionDelay,
                    policy,
                    bailInCoveragePercent,
                    bailOutCoveragePercent,
                    localBankLiquidity,
                    eLocalBankDebts,
                    bankDistressCount,
                    bailMoneyThisStep,
                    deficitThisStep
                );

                // F.3: Update total bail money
                // After the bail money for this timestep is finalized,
                // update the total tracker
                if (policy == 0) {
                    totalBailInMoney += bailMoneyThisStep;
                }
                else {
                    totalBailOutMoney += bailMoneyThisStep;
                }

                // F.6: Update total grace money
                // After all banks have been processed, add the grace money
                // used this timestep to the total grace money tracker
                totalDeficit += deficitThisStep;


                if (step == timesteps - 1) {
                    uint64_t reducedDeficitMoney = 0;

                    MPI_Reduce(
                        &totalDeficit,
                        &reducedDeficitMoney,
                        1,
                        MPI_UNSIGNED_LONG_LONG,
                        MPI_SUM,
                        0,
                        MPI_COMM_WORLD
                    );
                }
            }
            auto timestepEnd = chrono::high_resolution_clock::now();
            uint64_t elapsedMsTimestep = chrono::duration_cast<chrono::milliseconds>(
                timestepEnd - timestepStart
            ).count();

            if (mpiRank == 0) {
                printf(
                    "elapsedMs=%llu (TIMESTEP LOOP)\n",
                    (unsigned long long)elapsedMsTimestep
                );
                fflush(
                    stdout
                );
            }

            // Build final state components after timed loop ends
            long long localBankLiquiditySum = 0;
            long long localFirmLiquiditySum = 0;
            unsigned long long localBankDebtSum = 0ULL;
            unsigned long long localFirmDebtSum = 0ULL;
            unsigned long long localBankDistressSum = 0ULL;
            unsigned long long localBailMoneySum = 0ULL;
            unsigned long long localDeficitSum = 0ULL;

            for (unsigned int bankIndex = 0;
                bankIndex < bankCountForRank;
                bankIndex++) {

                localBankLiquiditySum +=
                    (long long)localBankLiquidity[bankIndex];

                localBankDistressSum +=
                    (unsigned long long)bankDistressCount[bankIndex];

                for (const DebtEntry& debt : eLocalBankDebts[bankIndex]) {
                    localBankDebtSum +=
                        (unsigned long long)debt.amount;
                }
            }

            for (unsigned int firmIndex = 0;
                firmIndex < firmCountForRank;
                firmIndex++) {

                localFirmLiquiditySum +=
                    (long long)localFirmLiquidity[firmIndex];

                for (const DebtEntry& debt : localFirmDebts[firmIndex]) {
                    localFirmDebtSum +=
                        (unsigned long long)debt.amount;
                }
            }

            if (policy == 0) {
                localBailMoneySum =
                    (unsigned long long)totalBailInMoney;
            }
            else {
                localBailMoneySum =
                    (unsigned long long)totalBailOutMoney;
            }

            localDeficitSum =
                (unsigned long long)totalDeficit;

            long long globalBankLiquiditySum = 0;
            long long globalFirmLiquiditySum = 0;
            unsigned long long globalBankDebtSum = 0ULL;
            unsigned long long globalFirmDebtSum = 0ULL;
            unsigned long long globalBankDistressSum = 0ULL;
            unsigned long long globalBailMoneySum = 0ULL;
            unsigned long long globalDeficitSum = 0ULL;

            MPI_Reduce(
                &localBankLiquiditySum,
                &globalBankLiquiditySum,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            MPI_Reduce(
                &localFirmLiquiditySum,
                &globalFirmLiquiditySum,
                1,
                MPI_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            MPI_Reduce(
                &localBankDebtSum,
                &globalBankDebtSum,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            MPI_Reduce(
                &localFirmDebtSum,
                &globalFirmDebtSum,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            MPI_Reduce(
                &localBankDistressSum,
                &globalBankDistressSum,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            MPI_Reduce(
                &localBailMoneySum,
                &globalBailMoneySum,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            MPI_Reduce(
                &localDeficitSum,
                &globalDeficitSum,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                0,
                MPI_COMM_WORLD
            );

            if (mpiRank == 0) {
                // Combine final state totals into one validation checksum
                auto mixChecksum = [](
                    uint64_t checksum,
                    uint64_t value
                ) {
                    checksum ^= value;
                    checksum *= 1099511628211ULL;
                    return checksum;
                };

                uint64_t finalChecksum = 1469598103934665603ULL;

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalBankLiquiditySum
                );

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalFirmLiquiditySum
                );

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalBankDebtSum
                );

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalFirmDebtSum
                );

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalBankDistressSum
                );

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalBailMoneySum
                );

                finalChecksum = mixChecksum(
                    finalChecksum,
                    (uint64_t)globalDeficitSum
                );

                printf(
                    "checksum run=%d policy=%d "
                    "bankLiquidity=%lld firmLiquidity=%lld "
                    "bankDebt=%llu firmDebt=%llu "
                    "bankDistress=%llu bailMoney=%llu "
                    "deficit=%llu finalChecksum=%llu\n",
                    run,
                    policy,
                    globalBankLiquiditySum,
                    globalFirmLiquiditySum,
                    globalBankDebtSum,
                    globalFirmDebtSum,
                    globalBankDistressSum,
                    globalBailMoneySum,
                    globalDeficitSum,
                    (unsigned long long)finalChecksum
                );

                fflush(
                    stdout
                );
            }

        }
    }

}
