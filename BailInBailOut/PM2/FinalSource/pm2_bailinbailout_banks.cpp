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

 // Build an MPI datatype for fixed layout message structs
MPI_Datatype BailInBailOut::makeStructType(
    int count,
    const MPI_Aint* offsets,
    const MPI_Datatype* types,
    MPI_Aint extent
) {
    vector<int> blockLengths(
        count,
        1
    );
    MPI_Datatype raw;
    MPI_Type_create_struct(
        count,
        blockLengths.data(),
        offsets,
        types,
        &raw
    );
    MPI_Datatype resized;
    MPI_Type_create_resized(
        raw,
        0,
        extent,
        &resized
    );
    MPI_Type_commit(
        &resized
    );
    MPI_Type_free(
        &raw
    );
    return resized;
}

// Create datatype for bank liquidity delta messages
MPI_Datatype BailInBailOut::c1MakeBankDeltaMessageType() {
    MPI_Aint offsets[3] = {
        (MPI_Aint)offsetof(
            c1BankDeltaMessage,
            bankGlobalId
        ),
        (MPI_Aint)offsetof(
            c1BankDeltaMessage,
            pad
        ),
        (MPI_Aint)offsetof(
            c1BankDeltaMessage,
            delta
        )
    };
    MPI_Datatype types[3] = {
        MPI_UNSIGNED,
        MPI_UNSIGNED,
        MPI_LONG_LONG
    };
    return makeStructType(
        3,
        offsets,
        types,
        static_cast<MPI_Aint>(
            sizeof(
                c1BankDeltaMessage
                )
            )
    );
}

// Route wage deposits and firm repayments to owning bank ranks
void BailInBailOut::c1SendMoneyToBanks(
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
) {
    // allToAllvExchange clears c1ToRank inner vectors at the end of its last
    // call, so inner capacity persists across timesteps without reallocation.
    // Send wage deposits for banks that received worker income
    for (uint32_t bankGID : wageNonZeroIds) {
        int64_t delta = wageDepositBuffer[bankGID] + bankIncomingFromFirmRepay[bankGID];
        unsigned int ownerRank = ownerRankFromGlobalId(
            bankCountTotal,
            mpiSize,
            bankGID
        );
        c1ToRank[ownerRank].push_back(
            c1BankDeltaMessage{ bankGID, 0, delta }
        );
    }
    // Send firm repayments that were not already covered by wages
    for (uint32_t bankGID : dirtyBankRepayIds) {
        if (wageDepositBuffer[bankGID] != 0) {
            continue;
        }
        int64_t delta = bankIncomingFromFirmRepay[bankGID];
        if (delta == 0) {
            continue;
        }
        unsigned int ownerRank = ownerRankFromGlobalId(
            bankCountTotal,
            mpiSize,
            bankGID
        );
        c1ToRank[ownerRank].push_back(
            c1BankDeltaMessage{ bankGID, 0, delta }
        );
    }
    for (uint32_t id : dirtyBankRepayIds) {
        bankIncomingFromFirmRepay[id] = 0;
    }
    dirtyBankRepayIds.clear();
    // Exchange all bank deltas and apply only messages owned here
    static MPI_Datatype bankDeltaType = c1MakeBankDeltaMessageType();
    static vector<c1BankDeltaMessage> recvBuffer;
    allToAllvExchange(
        mpiRank,
        mpiSize,
        c1ToRank,
        recvBuffer,
        bankDeltaType
    );
    for (unsigned int i = 0; i < recvBuffer.size(); i++) {
        unsigned int localIndex = recvBuffer[i].bankGlobalId - bankGlobalStartIndex;
        localBankLiquidity[localIndex] += recvBuffer[i].delta;
    }
}

// Create datatype for firm loan request messages
MPI_Datatype BailInBailOut::d1MakeFirmLoanRequestType() {
    MPI_Aint offsets[6] = {
        (MPI_Aint)offsetof(
            FirmLoanRequest,
            firmGlobalId
        ),
        (MPI_Aint)offsetof(
            FirmLoanRequest,
            lenderBankGlobalId
        ),
        (MPI_Aint)offsetof(
            FirmLoanRequest,
            amountRequested
        ),
        (MPI_Aint)offsetof(
            FirmLoanRequest,
            lenderInterestRate
        ),
        (MPI_Aint)offsetof(
            FirmLoanRequest,
            pad16
        ),
        (MPI_Aint)offsetof(
            FirmLoanRequest,
            pad32
        )
    };
    MPI_Datatype types[6] = {
        MPI_UNSIGNED,
        MPI_UNSIGNED,
        MPI_UNSIGNED_LONG_LONG,
        MPI_UNSIGNED_SHORT,
        MPI_UNSIGNED_SHORT,
        MPI_UNSIGNED
    };
    return makeStructType(
        6,
        offsets,
        types,
        static_cast<MPI_Aint>(
            sizeof(
                FirmLoanRequest
                )
            )
    );
}


// Route firm loan requests to lender bank owner ranks
void BailInBailOut::d1SendFirmLoanRequests(
    unsigned int mpiRank,
    unsigned int mpiSize,
    vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank,
    vector<FirmLoanRequest>& receivedLoanRequests
) {
    static MPI_Datatype firmLoanReqType = d1MakeFirmLoanRequestType();
    allToAllvExchange(
        mpiRank,
        mpiSize,
        firmLoanRequestsToRank,
        receivedLoanRequests,
        firmLoanReqType
    );
}



// Process incoming firm loan requests on lender bank owners
void BailInBailOut::d2ProcessFirmLoanRequests(
    unsigned int firmCountTotal,
    unsigned int bankCountForRank,
    unsigned int bankGlobalStartIndex,
    unsigned int mpiSize,
    unsigned short maxFirmLoanPercent,
    vector<int64_t>& localBankLiquidity,
    vector<FirmLoanRequest>& receivedLoanRequests,
    vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
    vector<uint64_t>& bankFirmLoanOutflow
) {
    for (unsigned int r = 0; r < mpiSize; r++) {
        firmLoanAcceptancesToRank[r].clear();
    }
    fill(
        bankFirmLoanOutflow.begin(),
        bankFirmLoanOutflow.end(),
        0
    );

    // Reuse per-bank request buckets across timesteps
    static vector<vector<FirmLoanRequest>> requestsByLocalBank;
    if (requestsByLocalBank.size() != bankCountForRank) {
        requestsByLocalBank.clear();
        requestsByLocalBank.resize(
            bankCountForRank
        );
    }
    else {
        for (unsigned int b = 0; b < bankCountForRank; b++) {
            requestsByLocalBank[b].clear();
        }
    }

    // Group requests by the local lender bank that owns them
    for (unsigned int i = 0; i < receivedLoanRequests.size(); i++) {
        FirmLoanRequest& req = receivedLoanRequests[i];

        unsigned int localBankIndex = req.lenderBankGlobalId - bankGlobalStartIndex;

        if (localBankIndex < bankCountForRank) {
            requestsByLocalBank[localBankIndex].push_back(
                req
            );
        }
    }

    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {
        vector<FirmLoanRequest>& reqs = requestsByLocalBank[localBank];
        if (reqs.empty()) {
            continue;
        }

        // Stable deterministic order before proportional granting
        sort(reqs.begin(), reqs.end(),
            [](
                const FirmLoanRequest& a,
                const FirmLoanRequest& b
                ) {
                    if (a.firmGlobalId != b.firmGlobalId) {
                        return a.firmGlobalId < b.firmGlobalId;
                    }
                    if (a.amountRequested != b.amountRequested) {
                        return a.amountRequested < b.amountRequested;
                    }
                    return a.lenderBankGlobalId < b.lenderBankGlobalId;
            }
        );

        vector<uint64_t> grants(
            reqs.size(),
            0
        );

        int64_t liqSigned = localBankLiquidity[localBank];
        if (liqSigned > 0) {
            uint64_t liquidity = (uint64_t)liqSigned;
            uint64_t stepCap = (liquidity * (uint64_t)maxFirmLoanPercent) / 100ULL;
            uint64_t remainingToLend = (stepCap < liquidity) ? stepCap : liquidity;

            if (remainingToLend > 0) {
                uint64_t totalRequested = 0;
                for (size_t k = 0; k < reqs.size(); k++) {
                    totalRequested += reqs[k].amountRequested;
                }

                if (totalRequested > 0) {
                    vector<uint64_t> remainders(
                        reqs.size(),
                        0
                    );
                    uint64_t sumGrants = 0;

                    // Allocate lender capacity proportionally by request size
                    for (size_t k = 0; k < reqs.size(); k++) {
                        unsigned __int128 scaled =
                            (unsigned __int128)remainingToLend * reqs[k].amountRequested;

                        grants[k] = (uint64_t)(scaled / totalRequested);
                        remainders[k] = (uint64_t)(scaled % totalRequested);

                        if (grants[k] > reqs[k].amountRequested) {
                            grants[k] = reqs[k].amountRequested;
                        }

                        sumGrants += grants[k];
                    }

                    // Floor division may leave units to distribute
                    uint64_t leftover = remainingToLend - sumGrants;
                    if (leftover > 0) {
                        vector<size_t> order;
                        order.reserve(
                            reqs.size()
                        );

                        for (size_t k = 0; k < reqs.size(); k++) {
                            if (reqs[k].amountRequested > grants[k]) {
                                order.push_back(
                                    k
                                );
                            }
                        }

                        sort(order.begin(), order.end(),
                            [&](
                                size_t a,
                                size_t b
                                ) {
                                    if (remainders[a] != remainders[b]) {
                                        return remainders[a] > remainders[b];
                                    }
                                    return reqs[a].firmGlobalId < reqs[b].firmGlobalId;
                            }
                        );

                        for (
                            size_t i = 0;
                            i < order.size() && leftover > 0;
                            i++
                            ) {
                            grants[order[i]]++;
                            leftover--;
                        }
                    }
                }
            }
        }

        // Return acceptance messages to the borrowing firm owners
        for (size_t k = 0; k < reqs.size(); k++) {
            unsigned int firmOwner = ownerRankFromGlobalId(
                firmCountTotal,
                mpiSize,
                reqs[k].firmGlobalId
            );

            firmLoanAcceptancesToRank[firmOwner].push_back(
                d2FirmLoanAcceptance{
                    reqs[k].firmGlobalId,
                    reqs[k].lenderBankGlobalId,
                    grants[k]
                }
            );

            if (grants[k] != 0) {
                bankFirmLoanOutflow[localBank] += grants[k];
            }
        }
    }

    receivedLoanRequests.clear();
}

// Create datatype for firm loan acceptance messages
MPI_Datatype BailInBailOut::d4MakeFirmLoanAcceptanceType() {
    MPI_Aint offsets[3] = {
        (MPI_Aint)offsetof(
            d2FirmLoanAcceptance,
            firmGlobalId
        ),
        (MPI_Aint)offsetof(
            d2FirmLoanAcceptance,
            lenderBankGlobalId
        ),
        (MPI_Aint)offsetof(
            d2FirmLoanAcceptance,
            amountGranted
        )
    };
    MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };
    return makeStructType(
        3,
        offsets,
        types,
        static_cast<MPI_Aint>(
            sizeof(
                d2FirmLoanAcceptance
                )
            )
    );
}

// Route firm loan acceptances back to borrower firm ranks
void BailInBailOut::d4SendFirmLoanAcceptances(
    unsigned int mpiRank,
    unsigned int mpiSize,
    vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
    vector<d2FirmLoanAcceptance>& receivedAcceptances
) {
    static MPI_Datatype firmLoanAccType = d4MakeFirmLoanAcceptanceType();
    allToAllvExchange(
        mpiRank,
        mpiSize,
        firmLoanAcceptancesToRank,
        receivedAcceptances,
        firmLoanAccType
    );
}

// Apply granted firm loans and record borrower side debt
void BailInBailOut::d5ApplyFirmLoanAcceptances(
    unsigned int firmGlobalStartIndex,
    vector<int64_t>& localFirmLiquidity,
    vector<vector<DebtEntry>>& localFirmDebts,
    vector<d2FirmLoanAcceptance>& receivedAcceptances
) {
    for (unsigned int i = 0; i < receivedAcceptances.size(); i++) {

        d2FirmLoanAcceptance& acc = receivedAcceptances[i];

        // Convert global firm id to local firm index
        unsigned int localFirmIndex = acc.firmGlobalId - firmGlobalStartIndex;

        // If routing is correct, this should always be in range
        if (localFirmIndex >= localFirmLiquidity.size()) {
            continue;
        }

        if (acc.amountGranted == 0) {
            continue;
        }

        // Apply cash inflow to firm
        localFirmLiquidity[localFirmIndex] += (int64_t)acc.amountGranted;

        // Record debt on the borrower side
        localFirmDebts[localFirmIndex].push_back(
            DebtEntry{
                acc.lenderBankGlobalId,
                acc.amountGranted
            }
        );
    }

    // Inbox is one step use
    receivedAcceptances.clear();
}

// Create datatype for interbank repayment messages
MPI_Datatype BailInBailOut::e1MakeInterbankRepaymentType() {
    MPI_Aint offsets[3] = {
        (MPI_Aint)offsetof(
            e1InterbankRepaymentMessage,
            lenderBankGlobalId
        ),
        (MPI_Aint)offsetof(
            e1InterbankRepaymentMessage,
            pad
        ),
        (MPI_Aint)offsetof(
            e1InterbankRepaymentMessage,
            amount
        )
    };
    MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };
    return makeStructType(
        3,
        offsets,
        types,
        static_cast<MPI_Aint>(
            sizeof(
                e1InterbankRepaymentMessage
                )
            )
    );
}

// Repay interbank debt from borrower banks to lender owners
void BailInBailOut::e1RepayInterbankLoans(
    unsigned int bankCountTotal,
    unsigned int mpiRank,
    unsigned int mpiSize,
    unsigned int bankGlobalStartIndex,
    unsigned int bankCountForRank,
    unsigned short bankRepayPercent,
    vector<int64_t>& localBankLiquidity,
    vector<vector<DebtEntry>>& localBankDebts,
    vector<vector<e1InterbankRepaymentMessage>>& e1ToRank
) {
    // e1ToRank is caller-owned; allToAllvExchange clears its inner vectors
    // at the end so push_back here reuses the persistent capacity.
    vector<vector<e1InterbankRepaymentMessage>>& toRank = e1ToRank;

    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {

        // Skip repayment if borrower bank has no positive liquidity
        if (localBankLiquidity[localBank] <= 0) {
            continue;
        }

        vector<DebtEntry>& debts = localBankDebts[localBank];

        // Skip if nothing to repay
        if (debts.empty()) {
            continue;
        }

        for (unsigned int i = 0; i < debts.size(); i++) {
            DebtEntry& entry = debts[i];

            // Skip this debt entry if it's 0
            if (entry.amount == 0) {
                continue;
            }

            // bankRepayPercent% of remaining debt
            uint64_t repay =
                (uint64_t)((entry.amount * (uint64_t)bankRepayPercent) / 100ULL);

            if (repay == 0) {
                continue;
            }

            // Don't repay more cash than borrower has
            uint64_t cash =
                (localBankLiquidity[localBank] > 0) ?
                (uint64_t)localBankLiquidity[localBank] : 0ULL;

            // If there's not enough to pay off the desired amount,
            // pay all of the available cash
            if (repay > cash) {
                repay = cash;
            }

            // If there's nothing left to repay, stop
            if (repay == 0) {
                break;
            }

            // Borrower subtracts instantly
            localBankLiquidity[localBank] -= (int64_t)repay;

            // Reduce the borrower's stored liability
            entry.amount -= (uint64_t)repay;

            // Send credit to lender owner rank similarly to C.1
            unsigned int lenderGlobalId = entry.lenderBankGlobalId;

            unsigned int ownerRank = ownerRankFromGlobalId(
                bankCountTotal,
                mpiSize,
                lenderGlobalId
            );

            // Route repayment credit to lender bank owner
            toRank[ownerRank].push_back(e1InterbankRepaymentMessage{
                lenderGlobalId,
                0,
                repay
                });

            // If borrower cash is now below 0, stop repaying further debts
            if (localBankLiquidity[localBank] <= 0) {
                break;
            }
        }

        // Remove fully repaid entries to save space
        unsigned int write = 0;
        for (unsigned int i = 0; i < debts.size(); i++) {
            if (debts[i].amount != 0) {
                debts[write++] = debts[i];
            }
        }
        debts.resize(
            write
        );
    }

    // Exchange and apply credits on lender owners
    static MPI_Datatype msgType = e1MakeInterbankRepaymentType();



    static vector<e1InterbankRepaymentMessage> received;
    allToAllvExchange<e1InterbankRepaymentMessage>(
        mpiRank,
        mpiSize,
        toRank,
        received,
        msgType
    );

    // Apply lender credits to banks owned by this rank
    for (unsigned int i = 0; i < received.size(); i++) {
        unsigned int lenderGlobalId = received[i].lenderBankGlobalId;
        uint64_t amount = received[i].amount;

        unsigned int lenderLocal = lenderGlobalId - bankGlobalStartIndex;
        if (lenderLocal < bankCountForRank) {
            localBankLiquidity[lenderLocal] += (int64_t)amount;
        }
        else {
            // Should NEVER be entered. Indicates mismatch in range mapping.
            printf(
                "E1 ERROR: Received lenderGlobalId=%u not owned by this rank\n",
                lenderGlobalId
            );
        }
    }
}

// Apply interest growth to both bank and firm debt books
void BailInBailOut::e2ApplyInterestOnAllLoans(
    vector<vector<DebtEntry>>& eLocalBankDebts,
    vector<vector<DebtEntry>>& localFirmDebts,
    uint64_t baseSeed,
    unsigned int run,
    uint64_t minVal,
    uint64_t maxVal
) {

    // Go through every BANK debt entry
    for (int i = 0; i < eLocalBankDebts.size(); i++) {
        for (DebtEntry& entry : eLocalBankDebts[i]) {
            entry.amount = (entry.amount * (100 +
                getInterestRate(
                    baseSeed,
                    run,
                    entry.lenderBankGlobalId,
                    minVal,
                    maxVal
                ))) / 100ULL;
        }
    }

    // Go through every FIRM debt entry
    for (int i = 0; i < localFirmDebts.size(); i++) {
        for (DebtEntry& entry : localFirmDebts[i]) {
            entry.amount = (entry.amount * (100 +
                getInterestRate(
                    baseSeed,
                    run,
                    entry.lenderBankGlobalId,
                    minVal,
                    maxVal
                ))) / 100ULL;
        }
    }

}

// Create datatype for interbank loan request messages
MPI_Datatype BailInBailOut::e3MakeInterbankLoanRequestType() {
    MPI_Aint offsets[6] = {
        (MPI_Aint)offsetof(
            e3InterbankLoanRequest,
            borrowerBankGlobalId
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanRequest,
            lenderBankGlobalId
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanRequest,
            amountRequested
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanRequest,
            lenderInterestRate
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanRequest,
            pad16
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanRequest,
            pad32
        )
    };
    MPI_Datatype types[6] = {
        MPI_UNSIGNED,
        MPI_UNSIGNED,
        MPI_UNSIGNED_LONG_LONG,
        MPI_UNSIGNED_SHORT,
        MPI_UNSIGNED_SHORT,
        MPI_UNSIGNED
    };
    return makeStructType(
        6,
        offsets,
        types,
        static_cast<MPI_Aint>(
            sizeof(
                e3InterbankLoanRequest
                )
            )
    );
}

// Create datatype for interbank loan acceptance messages
MPI_Datatype BailInBailOut::e3MakeInterbankLoanAcceptanceType() {
    MPI_Aint offsets[3] = {
        (MPI_Aint)offsetof(
            e3InterbankLoanAcceptance,
            borrowerBankGlobalId
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanAcceptance,
            lenderBankGlobalId
        ),
        (MPI_Aint)offsetof(
            e3InterbankLoanAcceptance,
            amountGranted
        )
    };
    MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };
    return makeStructType(
        3,
        offsets,
        types,
        static_cast<MPI_Aint>(
            sizeof(
                e3InterbankLoanAcceptance
                )
            )
    );
}







// Wrapper keeps phase E.3 call site independent of internal strategy
void BailInBailOut::e3BorrowInterbankIfNegativeLiquidity(
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
) {
    e3BorrowInterbankIfNegativeLiquidityDense(
        bankCountTotal,
        mpiRank,
        mpiSize,
        bankGlobalStartIndex,
        bankCountForRank,
        step,
        maxInterbankLenderSamplingK,
        maxInterbankLoanPercent,
        baseSeed,
        run,
        minInterestRate,
        maxInterestRate,
        localBankNeighbors,
        bankAllgatherRecvCounts,
        bankAllgatherDispls,
        localBankLiquidity,
        localBankDebts,
        lenderLiquidityView,
        samplingIndices,
        reqsByLocalLender,
        e3RequestsToRank,
        e3AcceptancesToRank
    );
}

// Dense E.3 path samples lenders using a full liquidity view
void BailInBailOut::e3BorrowInterbankIfNegativeLiquidityDense(
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
) {
    // E.3.1 gathers current bank liquidity for lender capacity estimates
    //
    // Every rank needs the current liquidity of every bank that could be a
    // lender. With interbank density ~21% and 71,680 banks, essentially every
    // rank's neighbor set touches every other rank, so the prior query based
    // alltoallv was functionally an all-gather with O(mpiSize^2) message
    // startup cost. MPI_Allgatherv delivers the same data in one collective
    // with O(log mpiSize) latency.
    //
    // bankAllgatherRecvCounts[r] and bankAllgatherDispls[r] describe how
    // computeRange partitions bankCountTotal across ranks (supports uneven
    // splits where the first `remainder` ranks carry one extra bank).
    // lenderLiquidityView ends up in global-bank-ID order across all ranks,
    // so it indexes directly by global bank ID without any INT64_MIN sentinel.

    if (mpiSize > 1) {
        MPI_Allgatherv(
            localBankLiquidity.data(),
            (int)bankCountForRank,
            MPI_LONG_LONG,
            lenderLiquidityView.data(),
            bankAllgatherRecvCounts.data(),
            bankAllgatherDispls.data(),
            MPI_LONG_LONG,
            MPI_COMM_WORLD
        );
    }

    // E.3.2 builds interbank loan requests for negative liquidity banks

    vector<vector<e3InterbankLoanRequest>>& requestsToRank = e3RequestsToRank;

    // Convert observed lender liquidity into max lendable amount
    auto getLenderCapFromLiquidity = [&](
        int64_t lenderLiquidity
        ) -> uint64_t {
            if (lenderLiquidity <= 0) {
                return 0ULL;
            }
            uint64_t liq = (uint64_t)lenderLiquidity;
            return (
                liq * (uint64_t)maxInterbankLoanPercent
                ) / 100ULL;
        };

    // Reused across borrowers: stores at most maxInterbankLenderSamplingK
    // (i, j) pairs from the partial Fisher-Yates so we can undo them.
    static vector<pair<unsigned int, unsigned int>> samplingSwaps;
    samplingSwaps.clear();
    if (samplingSwaps.capacity() < (size_t)maxInterbankLenderSamplingK) {
        samplingSwaps.reserve(
            (size_t)maxInterbankLenderSamplingK
        );
    }

    // Only banks with negative liquidity try to borrow interbank funds
    for (unsigned int localBorrower = 0; localBorrower < bankCountForRank; localBorrower++) {
        // If liquidity is not negative, don't need loans
        if (localBankLiquidity[localBorrower] >= 0) {
            continue;
        }

        // Get the global ID of the borrower bank
        unsigned int borrowerGlobalId = bankGlobalStartIndex + localBorrower;

        // Figure out the deficit
        // Deficit is how much is needed to be loaned
        uint64_t deficitRemaining = (uint64_t)(-localBankLiquidity[localBorrower]);

        // Extra check just to make sure there's some deficit
        if (deficitRemaining == 0) {
            continue;
        }

        // All the candidates for loaner
        vector<unsigned int>& candidates = localBankNeighbors[localBorrower];

        // If this bank how no loaner candidates, skip
        if (candidates.empty()) {
            continue;
        }

        unsigned int sampleCount = (unsigned int)maxInterbankLenderSamplingK;

        // If too little nodes to reach the max samplers, then
        // set the sampling size to total candidate count
        if (sampleCount > (
            unsigned int
        )candidates.size()) {
            sampleCount = (unsigned int)candidates.size();
        }

        // Sampled lenders are ranked before requests are created
        vector<e3LenderOption> lenderOptions;
        lenderOptions.reserve(
            sampleCount
        );

        {
            // Partial Fisher-Yates over the pre-allocated identity-init
            // samplingIndices buffer. Record each (i, j) swap and reverse
            // them after sampling so the buffer returns to identity for
            // the next borrower.
            samplingSwaps.clear();
            Xoshiro256 xRng(makeSeed(
                baseSeed,
                run,
                step,
                borrowerGlobalId,
                STREAM_INTERBANK_LENDER_SAMPLING
            ));
            for (unsigned int i = 0; i < sampleCount; i++) {
                uint64_t j = xRng.nextInRange(
                    (uint64_t)i,
                    (uint64_t)(candidates.size() - 1)
                );
                unsigned int jj = (unsigned int)j;
                swap(
                    samplingIndices[i],
                    samplingIndices[jj]
                );
                samplingSwaps.push_back(make_pair(
                    i,
                    jj
                ));
                unsigned int lenderGlobalId = candidates[samplingIndices[i]];

                unsigned short rate = getInterestRate(
                    baseSeed,
                    run,
                    lenderGlobalId,
                    minInterestRate,
                    maxInterestRate
                );

                int64_t liquidity = (mpiSize == 1) ?
                    localBankLiquidity[lenderGlobalId] :
                    lenderLiquidityView[lenderGlobalId];

                // Store lender score fields used by the sort below
                lenderOptions.push_back(e3LenderOption{
                    lenderGlobalId,
                    rate,
                    liquidity
                    });
            }
            // Restore identity by undoing swaps in reverse order
            for (size_t t = samplingSwaps.size(); t > 0; t--) {
                swap(
                    samplingIndices[samplingSwaps[t - 1].first],
                    samplingIndices[samplingSwaps[t - 1].second]
                );
            }
        }

        // Sort lenders by:
        // Lowest interest rate
        // If tied, higher liquidity
        // If STILL tied, lowest global ID
        sort(lenderOptions.begin(), lenderOptions.end(),
            [](
                const e3LenderOption& a,
                const e3LenderOption& b
                ) {
                    if (a.interestRate != b.interestRate) {
                        return a.interestRate < b.interestRate;
                    }
                    if (a.liquidity != b.liquidity) {
                        return a.liquidity > b.liquidity;
                    }
                    return a.lenderGlobalId < b.lenderGlobalId;
            }
        );

        // Request loans from multiple lenders if needed
        // Borrow as much as possible from best lender
        // If not enough, go to next, repeat
        for (unsigned int i = 0; i < lenderOptions.size(); i++) {
            if (deficitRemaining == 0) {
                break;
            }

            unsigned int lenderGlobalId = lenderOptions[i].lenderGlobalId;

            // If lender liquidity is unknown, still possible to request,
            // but it's hard to cap properly
            // In that case, request the remaining deficit, and let lender
            // handle whether it can grant that request or not
            int64_t lenderLiquidity = lenderOptions[i].liquidity;

            uint64_t plannedRequest = deficitRemaining;

            // Make sure the lender has some liquidity
            if (lenderLiquidity > 0) {
                // Planned request is capped by estimated lender capacity
                uint64_t estimatedCap = getLenderCapFromLiquidity(
                    lenderLiquidity
                );

                if (plannedRequest > estimatedCap) {
                    plannedRequest = estimatedCap;
                }
            }

            // If estimated cap is 0 (or lender had negative liquidity),
            // skip requesting this lender
            if (plannedRequest == 0) {
                continue;
            }

            // Get the owner rank of the lender
            unsigned int lenderOwnerRank = ownerRankFromGlobalId(
                bankCountTotal,
                mpiSize,
                lenderGlobalId
            );

            // Add this to pool of requests for this rank
            // Request is sent to the rank that owns the lender bank
            requestsToRank[lenderOwnerRank].push_back(e3InterbankLoanRequest{
                borrowerGlobalId,
                lenderGlobalId,
                plannedRequest,
                lenderOptions[i].interestRate
                });

            // Assume attempt to cover this portion of the deficit using this lender
            // If lender cannot fully grant, borrower will remain negative and
            // try again next timestep
            deficitRemaining -= plannedRequest;
        }
    }

    // E.3.3 routes loan requests to lender owning ranks

    static MPI_Datatype reqType = e3MakeInterbankLoanRequestType();

    static vector<e3InterbankLoanRequest> receivedRequests;
    allToAllvExchange(
        mpiRank,
        mpiSize,
        requestsToRank,
        receivedRequests,
        reqType,
        3301
    );

    // Lenders process requests deterministically and grant up to capacity

    // Reuse the persistent reqsByLocalLender buffer; clear inner vectors
    // (capacity preserved, so re-population is allocation-free after warmup).
    for (unsigned int b = 0; b < bankCountForRank; b++) {
        reqsByLocalLender[b].clear();
    }

    // Bucket received requests by local lender bank
    for (unsigned int i = 0; i < receivedRequests.size(); i++) {
        e3InterbankLoanRequest& req = receivedRequests[i];

        unsigned int localLender = req.lenderBankGlobalId - bankGlobalStartIndex;
        if (localLender < bankCountForRank) {
            reqsByLocalLender[localLender].push_back(
                req
            );
        }
    }

    vector<vector<e3InterbankLoanAcceptance>>& acceptancesToRank = e3AcceptancesToRank;

    // Each lender grants in deterministic request order up to capacity
    for (unsigned int localLender = 0; localLender < bankCountForRank; localLender++) {
        vector<e3InterbankLoanRequest>& reqs = reqsByLocalLender[localLender];
        if (reqs.empty()) {
            continue;
        }

        sort(reqs.begin(), reqs.end(),
            [](
                const e3InterbankLoanRequest& a,
                const e3InterbankLoanRequest& b
                ) {
                    if (a.borrowerBankGlobalId != b.borrowerBankGlobalId) {
                        return a.borrowerBankGlobalId < b.borrowerBankGlobalId;
                    }
                    if (a.amountRequested != b.amountRequested) {
                        return a.amountRequested < b.amountRequested;
                    }
                    return a.lenderBankGlobalId < b.lenderBankGlobalId;
            }
        );

        vector<uint64_t> grants(
            reqs.size(),
            0
        );

        int64_t lenderLiquiditySigned = localBankLiquidity[localLender];
        if (lenderLiquiditySigned > 0) {
            uint64_t lenderLiquidity = (uint64_t)lenderLiquiditySigned;
            uint64_t lendCap =
                (lenderLiquidity * (uint64_t)maxInterbankLoanPercent) / 100ULL;
            uint64_t remainingToLend =
                (lendCap < lenderLiquidity) ? lendCap : lenderLiquidity;

            for (unsigned int i = 0; i < reqs.size(); i++) {
                if (remainingToLend == 0) {
                    break;
                }

                uint64_t grant = reqs[i].amountRequested;
                if (grant > remainingToLend) {
                    grant = remainingToLend;
                }

                grants[i] = grant;
                remainingToLend -= grant;
            }
        }

        for (unsigned int i = 0; i < reqs.size(); i++) {
            if (grants[i] == 0) {
                continue;
            }

            unsigned int borrowerOwnerRank = ownerRankFromGlobalId(
                bankCountTotal,
                mpiSize,
                reqs[i].borrowerBankGlobalId
            );

            acceptancesToRank[borrowerOwnerRank].push_back(
                e3InterbankLoanAcceptance{
                    reqs[i].borrowerBankGlobalId,
                    reqs[i].lenderBankGlobalId,
                    grants[i]
                }
            );

            localBankLiquidity[localLender] -= (int64_t)grants[i];
        }
    }

    // Send granted interbank loan acceptances back to borrowers

    static MPI_Datatype accType = e3MakeInterbankLoanAcceptanceType();

    static vector<e3InterbankLoanAcceptance> receivedAcceptances;
    allToAllvExchange(
        mpiRank,
        mpiSize,
        acceptancesToRank,
        receivedAcceptances,
        accType,
        3501
    );

    // Borrowers apply granted loans and record interbank debt

    for (unsigned int i = 0; i < receivedAcceptances.size(); i++) {
        e3InterbankLoanAcceptance& acc = receivedAcceptances[i];

        unsigned int localBorrower = acc.borrowerBankGlobalId - bankGlobalStartIndex;
        if (localBorrower >= bankCountForRank) {
            continue;
        }

        if (acc.amountGranted == 0) {
            continue;
        }

        localBankLiquidity[localBorrower] += (int64_t)acc.amountGranted;

        localBankDebts[localBorrower].push_back(DebtEntry{
            acc.lenderBankGlobalId,
            acc.amountGranted
            });
    }
}


// Update distress counters and apply bail-in or bail-out policy
void BailInBailOut::f1f2UpdateBankDistressAndApplyIntervention(
    unsigned int bankCountForRank,
    unsigned int interventionDelay,
    unsigned int policy,
    unsigned short bailInCoveragePercent,
    unsigned short bailOutCoveragePercent,
    vector<int64_t>& localBankLiquidity,
    vector<vector<DebtEntry>>& localBankDebts,
    vector<unsigned int>& bankDistressCount,
    uint64_t& bailMoneyThisStep,
    uint64_t& deficitThisStep
) {
    bailMoneyThisStep = 0;
    deficitThisStep = 0;


    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {
        vector<DebtEntry>& debts = localBankDebts[localBank];
        // Compare liquidity against total interbank debt
        int64_t liq = localBankLiquidity[localBank];
        uint64_t totalDebtBefore = sumDebtU64(
            debts
        );

        bool isInsolvent;
        uint64_t deficitBefore;

        if (liq >= 0) {
            uint64_t liqU = (uint64_t)liq;
            isInsolvent = liqU < totalDebtBefore;
            deficitBefore = isInsolvent ? (totalDebtBefore - liqU) : 0;
        }
        else {
            isInsolvent = true;
            deficitBefore = totalDebtBefore + (uint64_t)(-liq);
        }

        // F.1: Update distress counter
        if (isInsolvent) {
            bankDistressCount[localBank]++;
        }
        // Reset if it's solvent
        else {
            bankDistressCount[localBank] = 0;
            continue;
        }

        // F.2: Apply intervention every interventionDelay distress timesteps
        unsigned int distressCount = bankDistressCount[localBank];

        bool shouldIntervene =
            (distressCount >= interventionDelay) &&
            (distressCount % interventionDelay == 0);
        // If it's not time to intervene yet, don't proceed
        if (!shouldIntervene) {
            continue;
        }

        // Policy 0 -> Bail-In
        if (policy == 0) {
            // Check if there is some debt
            if (totalDebtBefore > 0) {
                uint64_t targetDebtRemoval =
                    percentFloorU64(
                        deficitBefore,
                        bailInCoveragePercent
                    );

                // If the target removal is greater than what is needed,
                // only cover the total debt
                if (targetDebtRemoval > totalDebtBefore) {
                    targetDebtRemoval = totalDebtBefore;
                }

                if (targetDebtRemoval > 0) {
                    const size_t debtCount = debts.size();

                    // Amount of reduction applied to each debt entry
                    vector<uint64_t> debtReduction(
                        debtCount,
                        0
                    );
                    // Remainder used to distribute leftover reductions
                    vector<uint64_t> remainderByDebt(
                        debtCount,
                        0
                    );
                    // Total reduction assigned during the allocation pass
                    uint64_t assignedReduction = 0;

                    // Compute each debt reduction by liability share
                    for (size_t i = 0; i < debtCount; i++) {

                        // Skip debts that are already zero
                        if (debts[i].amount == 0) {
                            continue;
                        }

                        // Use 128 bit to avoid overflow when multiplying
                        unsigned __int128 scaled =
                            (unsigned __int128)targetDebtRemoval *
                            (unsigned __int128)debts[i].amount;

                        // Determine the proportional reduction
                        uint64_t reduction =
                            (uint64_t)(scaled / totalDebtBefore);
                        // Store remainder for distribution
                        uint64_t remainder =
                            (uint64_t)(scaled % totalDebtBefore);

                        // Never remove more than the actual debt amount
                        if (reduction > debts[i].amount) {
                            reduction = debts[i].amount;
                        }

                        debtReduction[i] = reduction;
                        remainderByDebt[i] = remainder;
                        assignedReduction += reduction;
                    }

                    // Remaining reduction still to be assigned
                    uint64_t leftover = targetDebtRemoval - assignedReduction;

                    if (leftover > 0) {
                        // Order of debts eligible for receiving the leftovers
                        vector<size_t> order;
                        order.reserve(
                            debtCount
                        );
                        // Only debts that still have remaining balance are eligible
                        for (size_t i = 0; i < debtCount; i++) {
                            if (debts[i].amount > debtReduction[i]) {
                                order.push_back(
                                    i
                                );
                            }
                        }

                        // Sort by remainder (largest first)
                        sort(
                            order.begin(),
                            order.end(),
                            [&](
                                size_t a,
                                size_t b
                                ) {
                                    if (remainderByDebt[a] != remainderByDebt[b]) {
                                        return remainderByDebt[a] > remainderByDebt[b];
                                    }
                                    if (debts[a].lenderBankGlobalId !=
                                        debts[b].lenderBankGlobalId) {
                                        return debts[a].lenderBankGlobalId <
                                            debts[b].lenderBankGlobalId;
                                    }
                                    return a < b;
                            }
                        );

                        // Distribute leftover reduction one unit at a time
                        for (
                            size_t k = 0;
                            k < order.size() && leftover > 0;
                            k++
                            ) {
                            size_t i = order[k];

                            if (debts[i].amount > debtReduction[i]) {
                                debtReduction[i] += 1;
                                leftover -= 1;
                            }
                        }
                    }

                    // Apply the reductions to the actual debt entries
                    uint64_t actuallyRemoved = 0;

                    for (size_t i = 0; i < debtCount; i++) {
                        uint64_t reduction = debtReduction[i];

                        if (reduction == 0) {
                            continue;
                        }

                        debts[i].amount -= reduction;
                        actuallyRemoved += reduction;
                    }

                    // Remove any debts that were fully eliminated
                    for (size_t i = 0; i < debts.size(); ) {
                        if (debts[i].amount == 0) {
                            debts[i] = debts.back();
                            debts.pop_back();
                        }
                        else {
                            i++;
                        }
                    }

                    // Track total debt removed by bail-in this timestep
                    bailMoneyThisStep += actuallyRemoved;
                }
            }
        }

        // Policy 1 -> Bail-Out
        else {
            if (bailOutCoveragePercent > 0) {
                // Find the injection amount as percentage of deficit
                uint64_t injection =
                    percentFloorU64(
                        deficitBefore,
                        bailOutCoveragePercent
                    );

                // Add that cash
                localBankLiquidity[localBank] += (int64_t)injection;
                // Track the liquidity
                bailMoneyThisStep += injection;
            }
        }

        // Track remaining deficit after intervention
        uint64_t totalDebtAfter = sumDebtU64(
            debts
        );

        bool stillInsolvent =
            ((long double)localBankLiquidity[localBank] < (long double)totalDebtAfter);

        if (stillInsolvent) {
            uint64_t deficitAfter =
                (uint64_t)((long double)totalDebtAfter -
                    (long double)localBankLiquidity[localBank]);

            deficitThisStep += deficitAfter;
        }
    }
}
