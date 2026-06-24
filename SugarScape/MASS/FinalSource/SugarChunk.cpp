#include "SugarChunk.h"
#include "MASS_base.h"

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>

// Map a global row to the chunk that owns it
static int rankForRow(
  int globalRow,
  int gridHeight,
  int numChunks
)
{
  if (numChunks <= 1) {
    return 0;
  }

  const int rowsPerRank = (gridHeight + numChunks - 1) / numChunks;
  int owner = globalRow / rowsPerRank;

  if (owner >= numChunks) {
    owner = numChunks - 1;
  }

  return owner;
}

// Get first global row owned by a chunk
static int firstRowForChunk(
  int chunk,
  int gridHeight,
  int numChunks
)
{
  const int rowsPerRank = (gridHeight + numChunks - 1) / numChunks;

  return chunk * rowsPerRank;
}

// Get one past last global row owned by a chunk
static int endRowForChunk(
  int chunk,
  int gridHeight,
  int numChunks
)
{
  const int rowsPerRank = (gridHeight + numChunks - 1) / numChunks;
  const int end = chunk * rowsPerRank + rowsPerRank;

  return std::min(
    end,
    gridHeight
  );
}

// Return movement direction priority for target tie breaking
static int directionIndexOf(
  int targetRow,
  int targetCol,
  int agentRow,
  int agentCol
)
{
  const int rowDiff = targetRow - agentRow;
  const int colDiff = targetCol - agentCol;

  if (rowDiff < 0 && colDiff == 0) {
    return 0;
  }

  if (rowDiff > 0 && colDiff == 0) {
    return 1;
  }

  if (rowDiff == 0 && colDiff > 0) {
    return 2;
  }

  if (rowDiff == 0 && colDiff < 0) {
    return 3;
  }

  return 4;
}

// Compute fixed sugar capacity for a global cell
static int sugarCapacityAt(
  int row,
  int col,
  int gridHeight,
  int gridWidth,
  int minCapacity,
  int maxCapacity
)
{
  const double peak1Row = gridHeight * 0.25;
  const double peak1Col = gridWidth * 0.25;
  const double peak2Row = gridHeight * 0.75;
  const double peak2Col = gridWidth * 0.75;
  const double maxRadius = std::min(
    static_cast<double>(
      gridHeight
      ),
    static_cast<double>(
      gridWidth
      )
  ) * 0.4;

  const double distance1 = std::sqrt(
    (row - peak1Row) * (row - peak1Row) +
    (col - peak1Col) * (col - peak1Col)
  );
  const double distance2 = std::sqrt(
    (row - peak2Row) * (row - peak2Row) +
    (col - peak2Col) * (col - peak2Col)
  );
  const double distance = std::min(
    distance1,
    distance2
  );

  if (distance <= maxRadius * 0.25) {
    return maxCapacity;
  }

  if (distance <= maxRadius * 0.50) {
    return std::min(
      maxCapacity,
      std::max(
        minCapacity,
        maxCapacity - 1
      )
    );
  }

  if (distance <= maxRadius * 0.75) {
    return std::min(
      maxCapacity,
      std::max(
        minCapacity,
        maxCapacity - 2
      )
    );
  }

  return minCapacity;
}

// Create a MASS Place instance
extern "C" Place* instantiate(
  void* argument
)
{
  return new SugarChunk(
    argument
  );
}

// Destroy a MASS Place instance
extern "C" void destroy(
  Place* object
)
{
  delete object;
}

SugarChunk::SugarChunk(
  void* argument
)
  : Place(
    argument
  ),
  chunkIndex_(
    -1
  ),
  numChunks_(
    1
  ),
  gridHeight_(
    0
  ),
  gridWidth_(
    0
  ),
  localStartRow_(
    0
  ),
  localRowCount_(
    0
  ),
  previousChunk_(
    -1
  ),
  nextChunk_(
    -1
  ),
  growthRate_(
    1
  ),
  sugarCapMin_(
    1
  ),
  sugarCapMax_(
    4
  ),
  metabolismMin_(
    1
  ),
  metabolismMax_(
    4
  ),
  visionMin_(
    1
  ),
  visionMax_(
    6
  ),
  initialWealth_(
    static_cast<unsigned int>(
      -1
      )
  ),
  outMessageBytes_(
    0
  ),
  currentStep_(
    0
  ),
  compoundStepPrimed_(
    false
  ),
  currentStamp_(
    1
  )
{
  outMessage = NULL;
  outMessage_size = 0;
}

// Phase tags stored at start of each boundary message
static const int PHASE_HALO_SUGAR = 1;
static const int PHASE_MOVE_CANDIDATES = 2;
static const int PHASE_WINNER_REPLIES = 3;
static const int PHASE_MIGRATIONS = 4;

// Each message starts with phase, above count, below count, padding
static const int HEADER_BYTES = 4 * sizeof(
  int
  );

// Allocate and clear shared boundary message buffer
void SugarChunk::ensureOutMessage()
{
  if (outMessage == NULL || outMessage_size != outMessageBytes_) {
    if (outMessage) {
      std::free(
        outMessage
      );
    }

    outMessage = std::malloc(
      outMessageBytes_
    );
    outMessage_size = outMessageBytes_;
  }

  std::memset(
    outMessage,
    0,
    outMessageBytes_
  );
}

// Reset message header for the next exchange phase
void SugarChunk::resetOutMessageHeader(
  int phase
)
{
  ensureOutMessage();

  int* header = static_cast<int*>(
    outMessage
    );
  header[0] = phase;
  header[1] = 0;
  header[2] = 0;
  header[3] = 0;
}

// Check whether this chunk owns a global row
bool SugarChunk::ownsRow(
  int globalRow
) const
{
  return globalRow >= localStartRow_ &&
    globalRow < localStartRow_ + localRowCount_;
}

// Find owner chunk for a global row
int SugarChunk::ownerChunkOfRow(
  int globalRow
) const
{
  return rankForRow(
    globalRow,
    gridHeight_,
    numChunks_
  );
}

// Return occupant id for a local cell
int SugarChunk::getOccupantAt(
  int globalRow,
  int col
) const
{
  if (
    !ownsRow(
      globalRow
    )
    ) {
    return -1;
  }

  return occupantIds_[
    (globalRow - localStartRow_) * gridWidth_ + col
  ];
}

// Set occupant id for a local cell
void SugarChunk::setOccupantAt(
  int globalRow,
  int col,
  int agentId
)
{
  if (
    !ownsRow(
      globalRow
    )
    ) {
    return;
  }

  occupantIds_[
    (globalRow - localStartRow_) * gridWidth_ + col
  ] = agentId;
}

// Clear occupant id for a local cell
void SugarChunk::clearOccupantAt(
  int globalRow,
  int col
)
{
  setOccupantAt(
    globalRow,
    col,
    -1
  );
}

// Rebuild local cell to agent lookup after migration
void SugarChunk::rebuildOccupantIndex()
{
  std::fill(
    occupantIds_.begin(),
    occupantIds_.end(),
    -1
  );

  for (const auto& agent : localAgents_) {
    if (
      ownsRow(
        agent.row
      ) &&
      agent.col >= 0 &&
      agent.col < gridWidth_
      ) {
      occupantIds_[
        (agent.row - localStartRow_) * gridWidth_ + agent.col
      ] = agent.id;
    }
  }
}

// Initialize chunk state, local sugar, and helper buffers
void* SugarChunk::init(
  void* argument
)
{
  if (!argument) {
    return NULL;
  }

  SugarChunkConfig* config = static_cast<SugarChunkConfig*>(
    argument
    );

  // Read benchmark configuration shared by all chunks
  chunkIndex_ = index[0];
  numChunks_ = size[0];
  gridHeight_ = config->height;
  gridWidth_ = config->width;
  growthRate_ = config->growthRate;
  sugarCapMin_ = config->sugarCapMin;
  sugarCapMax_ = config->sugarCapMax;
  metabolismMin_ = config->metabolismMin;
  metabolismMax_ = config->metabolismMax;
  visionMin_ = config->visionMin;
  visionMax_ = config->visionMax;
  initialWealth_ = config->initialWealth;
  outMessageBytes_ = config->outMessageBytes;
  reportAgentsCapacityBytes_ = config->reportAgentsBytes;
  reportSugarCapacityInts_ = config->reportSugarInts;

  // Prepare result buffers used by report calls
  reportAgentsBuf_.assign(
    reportAgentsCapacityBytes_,
    0
  );
  reportSugarBuf_.clear();

  // Compute owned row interval for this chunk
  localStartRow_ = firstRowForChunk(
    chunkIndex_,
    gridHeight_,
    numChunks_
  );

  const int localEndRow = endRowForChunk(
    chunkIndex_,
    gridHeight_,
    numChunks_
  );

  localRowCount_ = localEndRow - localStartRow_;

  // Store immediate neighbor ids for boundary exchange
  previousChunk_ = (
    chunkIndex_ > 0
    ? chunkIndex_ - 1
    : -1
    );
  nextChunk_ = (
    chunkIndex_ + 1 < numChunks_
    ? chunkIndex_ + 1
    : -1
    );

  // Allocate local grid data
  localSugar_.assign(
    localRowCount_ * gridWidth_,
    0
  );
  localCapacity_.assign(
    localRowCount_ * gridWidth_,
    0
  );
  sugarIsGrowing_.assign(
    static_cast<size_t>(
      localRowCount_
      ) * static_cast<size_t>(
        gridWidth_
        ),
    0
  );
  growingSugarCells_.clear();
  occupantIds_.assign(
    localRowCount_ * gridWidth_,
    -1
  );

  // Initialize sugar capacity and current sugar
  for (int row = 0; row < localRowCount_; ++row) {
    const int globalRow = localStartRow_ + row;

    for (int col = 0; col < gridWidth_; ++col) {
      const int capacity = sugarCapacityAt(
        globalRow,
        col,
        gridHeight_,
        gridWidth_,
        sugarCapMin_,
        sugarCapMax_
      );
      const int localIndex = row * gridWidth_ + col;
      int storedCapacity = capacity;

      if (storedCapacity < 0) {
        storedCapacity = 0;
      }

      if (storedCapacity > 255) {
        storedCapacity = 255;
      }

      localCapacity_[localIndex] = static_cast<uint8_t>(
        storedCapacity
        );
      localSugar_[localIndex] = storedCapacity;
    }
  }

  localAgents_.clear();
  ensureOutMessage();
  resetOutMessageHeader(
    PHASE_HALO_SUGAR
  );

  currentStep_ = 0;
  compoundStepPrimed_ = false;
  currentStamp_ = 1;

  const size_t localCellCount = static_cast<size_t>(
    localRowCount_
    ) * static_cast<size_t>(
      gridWidth_
      );

  // Allocate conflict resolution scratch buffers
  bestLocalAgentIndex_.assign(
    localCellCount,
    -1
  );
  bestAgentId_.assign(
    localCellCount,
    INT_MAX
  );
  cellStamp_.assign(
    localCellCount,
    0
  );

  return NULL;
}

// Keep only seed agents owned by this chunk
void* SugarChunk::seedInitialAgents(
  void* argument
)
{
  if (!argument) {
    return NULL;
  }

  SugarChunkSeedHeader* header = static_cast<SugarChunkSeedHeader*>(
    argument
    );
  AgentRec* agents = reinterpret_cast<AgentRec*>(
    static_cast<char*>(
      argument
      ) + sizeof(
        SugarChunkSeedHeader
        )
    );

  localAgents_.clear();
  localAgents_.reserve(
    header->numAgents
  );

  // Driver sends all agents, and each chunk filters its owned rows
  for (int i = 0; i < header->numAgents; ++i) {
    const AgentRec& agent = agents[i];

    if (
      ownsRow(
        agent.row
      ) &&
      agent.col >= 0 &&
      agent.col < gridWidth_
      ) {
      localAgents_.push_back(
        agent
      );
    }
  }

  rebuildOccupantIndex();

  return NULL;
}

// Seed redistribution is already handled by local filtering
void* SugarChunk::redistributeSeedAgents(
  void* /*argument*/
)
{
  return NULL;
}

// Mark a harvested cell so sparse regrowth visits it later
void SugarChunk::markSugarGrowing(
  int localCellIndex
)
{
  if (
    localCellIndex < 0 ||
    localCellIndex >= localRowCount_ * gridWidth_
    ) {
    return;
  }

  if (
    !sugarIsGrowing_[
      static_cast<size_t>(
        localCellIndex
        )
    ]
    ) {
    sugarIsGrowing_[
      static_cast<size_t>(
        localCellIndex
        )
    ] = 1;
    growingSugarCells_.push_back(
      localCellIndex
    );
  }
}

// Reset per cell conflict state lazily for current timestep
void SugarChunk::touchConflictCell(
  int localCellIndex
)
{
  if (
    cellStamp_[
      static_cast<size_t>(
        localCellIndex
        )
    ] != currentStamp_
    ) {
    cellStamp_[
      static_cast<size_t>(
        localCellIndex
        )
    ] = currentStamp_;
    bestLocalAgentIndex_[
      static_cast<size_t>(
        localCellIndex
        )
    ] = -1;
    bestAgentId_[
      static_cast<size_t>(
        localCellIndex
        )
    ] = INT_MAX;
  }
}

// Grow only cells that were harvested earlier
void* SugarChunk::growSugar(
  void* /*argument*/
)
{
  if (currentStep_ == 0 || growingSugarCells_.empty()) {
    return NULL;
  }

  size_t writeIndex = 0;

  // Compact list while removing cells that reached capacity
  for (
    size_t readIndex = 0;
    readIndex < growingSugarCells_.size();
    ++readIndex
    ) {
    const int localIndex = growingSugarCells_[readIndex];
    const int capacity = static_cast<int>(
      localCapacity_[localIndex]
      );
    const int grownSugar = localSugar_[localIndex] + growthRate_;

    if (grownSugar >= capacity) {
      localSugar_[localIndex] = capacity;
      sugarIsGrowing_[
        static_cast<size_t>(
          localIndex
          )
      ] = 0;
    }
    else {
      localSugar_[localIndex] = grownSugar;
      growingSugarCells_[writeIndex++] = localIndex;
    }
  }

  growingSugarCells_.resize(
    writeIndex
  );

  return NULL;
}

// Pack boundary sugar rows for neighboring chunks
void* SugarChunk::packHaloSugar(
  void* /*argument*/
)
{
  resetOutMessageHeader(
    PHASE_HALO_SUGAR
  );

  int topCount = std::min(
    visionMax_,
    localRowCount_
  );
  int bottomCount = std::min(
    visionMax_,
    localRowCount_
  );

  // Edge chunks don't send beyond grid boundary
  if (previousChunk_ < 0) {
    topCount = 0;
  }

  if (nextChunk_ < 0) {
    bottomCount = 0;
  }

  int* header = static_cast<int*>(
    outMessage
    );
  header[1] = topCount;
  header[2] = bottomCount;

  char* base = static_cast<char*>(
    outMessage
    ) + HEADER_BYTES;

  // Pack top owned rows for previous chunk
  if (topCount > 0) {
    std::memcpy(
      base,
      localSugar_.data(),
      topCount * gridWidth_ * sizeof(
        int
        )
    );
  }

  // Pack bottom owned rows for next chunk
  if (bottomCount > 0) {
    const int sourceOffset = (localRowCount_ - bottomCount) * gridWidth_;

    std::memcpy(
      base + topCount * gridWidth_ * sizeof(
        int
        ),
      localSugar_.data() + sourceOffset,
      bottomCount * gridWidth_ * sizeof(
        int
        )
    );
  }

  return NULL;
}

// Pick target cell for one agent using sugar, distance, then direction
void SugarChunk::computeTargetFor(
  const AgentRec& agent,
  const std::vector<int>& sugarBuffer,
  int sugarBufferStartRow,
  int sugarBufferRowCount,
  int& outRow,
  int& outCol
) const
{
  auto sugarAt = [&](
    int globalRow,
    int globalCol
    ) -> int {
      const int localRow = globalRow - sugarBufferStartRow;

      if (
        localRow < 0 ||
        localRow >= sugarBufferRowCount ||
        globalCol < 0 ||
        globalCol >= gridWidth_
        ) {
        return -1;
      }

      return sugarBuffer[localRow * gridWidth_ + globalCol];
    };

  int bestSugar = sugarAt(
    agent.row,
    agent.col
  );

  int candidateRows[100];
  int candidateCols[100];
  int candidateCount = 0;

  candidateRows[0] = agent.row;
  candidateCols[0] = agent.col;
  candidateCount = 1;

  // Scan visible cells in 4 cardinal directions
  for (int direction = 0; direction < 4; ++direction) {
    const int rowStep = SUGARCHUNK_DIRS[direction][0];
    const int colStep = SUGARCHUNK_DIRS[direction][1];

    for (int vision = 1; vision <= agent.vision; ++vision) {
      const int targetRow = agent.row + rowStep * vision;
      const int targetCol = agent.col + colStep * vision;

      if (
        targetRow < 0 ||
        targetRow >= gridHeight_ ||
        targetCol < 0 ||
        targetCol >= gridWidth_
        ) {
        break;
      }

      const int sugar = sugarAt(
        targetRow,
        targetCol
      );

      if (sugar < 0) {
        break;
      }

      // Track all cells with best visible sugar
      if (sugar > bestSugar) {
        bestSugar = sugar;
        candidateCount = 0;
        candidateRows[candidateCount] = targetRow;
        candidateCols[candidateCount] = targetCol;
        ++candidateCount;
      }
      else if (sugar == bestSugar && candidateCount < 100) {
        candidateRows[candidateCount] = targetRow;
        candidateCols[candidateCount] = targetCol;
        ++candidateCount;
      }
    }
  }

  if (candidateCount <= 1) {
    outRow = candidateRows[0];
    outCol = candidateCols[0];

    return;
  }

  int bestDistance = INT_MAX;
  int nearestRows[100];
  int nearestCols[100];
  int nearestCount = 0;

  // Keep only nearest cells among best sugar candidates
  for (int i = 0; i < candidateCount; ++i) {
    const int distance = std::abs(
      candidateRows[i] - agent.row
    ) + std::abs(
      candidateCols[i] - agent.col
    );

    if (distance < bestDistance) {
      bestDistance = distance;
      nearestCount = 0;
      nearestRows[nearestCount] = candidateRows[i];
      nearestCols[nearestCount] = candidateCols[i];
      ++nearestCount;
    }
    else if (distance == bestDistance && nearestCount < 100) {
      nearestRows[nearestCount] = candidateRows[i];
      nearestCols[nearestCount] = candidateCols[i];
      ++nearestCount;
    }
  }

  if (nearestCount <= 1) {
    outRow = nearestRows[0];
    outCol = nearestCols[0];

    return;
  }

  int bestIndex = 0;
  int bestDirection = directionIndexOf(
    nearestRows[0],
    nearestCols[0],
    agent.row,
    agent.col
  );

  // Break final ties by fixed direction order
  for (int i = 1; i < nearestCount; ++i) {
    const int direction = directionIndexOf(
      nearestRows[i],
      nearestCols[i],
      agent.row,
      agent.col
    );

    if (direction < bestDirection) {
      bestDirection = direction;
      bestIndex = i;
    }
  }

  outRow = nearestRows[bestIndex];
  outCol = nearestCols[bestIndex];
}

// Build local sugar view with neighbor halo rows
void* SugarChunk::computeTargets(
  void* /*argument*/
)
{
  const int haloCountAbove = std::min(
    visionMax_,
    localStartRow_
  );
  const int haloCountBelow = std::min(
    visionMax_,
    gridHeight_ - (localStartRow_ + localRowCount_)
  );

  const int sugarBufferStartRow = localStartRow_ - haloCountAbove;
  const int sugarBufferRowCount =
    haloCountAbove + localRowCount_ + haloCountBelow;

  std::vector<int> sugarBuffer(
    sugarBufferRowCount * gridWidth_,
    0
  );

  // Copy lower halo from previous chunk
  if (haloCountAbove > 0 && previousChunk_ >= 0) {
    int offset[1] = {
        -1
    };
    void* neighbor = getOutMessage(
      1,
      offset
    );

    if (neighbor) {
      int* header = static_cast<int*>(
        neighbor
        );
      const int neighborTopCount = header[1];
      const int neighborBottomCount = header[2];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES +
        neighborTopCount * gridWidth_ * sizeof(
          int
          );
      const int copyRows = std::min(
        haloCountAbove,
        neighborBottomCount
      );

      if (copyRows > 0) {
        const int sourceOffset =
          (neighborBottomCount - copyRows) * gridWidth_;

        std::memcpy(
          sugarBuffer.data() +
          (haloCountAbove - copyRows) * gridWidth_,
          base + sourceOffset * sizeof(
            int
            ),
          copyRows * gridWidth_ * sizeof(
            int
            )
        );
      }
    }
  }

  // Copy owned rows into center of sugar buffer
  std::memcpy(
    sugarBuffer.data() + haloCountAbove * gridWidth_,
    localSugar_.data(),
    localRowCount_ * gridWidth_ * sizeof(
      int
      )
  );

  // Copy upper halo from next chunk
  if (haloCountBelow > 0 && nextChunk_ >= 0) {
    int offset[1] = {
        1
    };
    void* neighbor = getOutMessage(
      1,
      offset
    );

    if (neighbor) {
      int* header = static_cast<int*>(
        neighbor
        );
      const int neighborTopCount = header[1];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES;
      const int copyRows = std::min(
        haloCountBelow,
        neighborTopCount
      );

      if (copyRows > 0) {
        std::memcpy(
          sugarBuffer.data() +
          (haloCountAbove + localRowCount_) * gridWidth_,
          base,
          copyRows * gridWidth_ * sizeof(
            int
            )
        );
      }
    }
  }

  const int localAgentCount = static_cast<int>(
    localAgents_.size()
    );

  targetCells_.assign(
    localAgentCount,
    {
        0,
        0
    }
  );
  moveWon_.assign(
    localAgentCount,
    0
  );

  // Compute requested target for every local agent
  for (int i = 0; i < localAgentCount; ++i) {
    int targetRow = 0;
    int targetCol = 0;

    computeTargetFor(
      localAgents_[i],
      sugarBuffer,
      sugarBufferStartRow,
      sugarBufferRowCount,
      targetRow,
      targetCol
    );

    targetCells_[i] = {
        targetRow,
        targetCol
    };
  }

  return NULL;
}

// Pack move requests whose target belongs to a neighbor chunk
void* SugarChunk::packMoveCandidates(
  void* /*argument*/
)
{
  sentCandidatesAbove_.clear();
  sentCandidatesBelow_.clear();
  sentCandidateIdxAbove_.clear();
  sentCandidateIdxBelow_.clear();

  const int localAgentCount = static_cast<int>(
    localAgents_.size()
    );

  for (int i = 0; i < localAgentCount; ++i) {
    const int targetRow = targetCells_[i].first;
    const int targetCol = targetCells_[i].second;

    if (
      ownsRow(
        targetRow
      )
      ) {
      continue;
    }

    MoveCandidateRec record{
        localAgents_[i].id,
        localAgents_[i].row,
        localAgents_[i].col,
        targetRow,
        targetCol
    };

    if (targetRow < localStartRow_) {
      sentCandidatesAbove_.push_back(
        record
      );
      sentCandidateIdxAbove_.push_back(
        i
      );
    }
    else {
      sentCandidatesBelow_.push_back(
        record
      );
      sentCandidateIdxBelow_.push_back(
        i
      );
    }
  }

  resetOutMessageHeader(
    PHASE_MOVE_CANDIDATES
  );

  int* header = static_cast<int*>(
    outMessage
    );
  const int sentUp = static_cast<int>(
    sentCandidatesAbove_.size()
    );
  const int sentDown = static_cast<int>(
    sentCandidatesBelow_.size()
    );

  header[1] = sentUp;
  header[2] = sentDown;

  const size_t recordBytes = sizeof(
    MoveCandidateRec
    );
  const size_t capacity = (outMessageBytes_ - HEADER_BYTES) / recordBytes;
  (void)capacity;

  char* base = static_cast<char*>(
    outMessage
    ) + HEADER_BYTES;

  // Above records are stored before below records
  if (sentUp > 0) {
    std::memcpy(
      base,
      sentCandidatesAbove_.data(),
      sentUp * recordBytes
    );
  }

  if (sentDown > 0) {
    std::memcpy(
      base + sentUp * recordBytes,
      sentCandidatesBelow_.data(),
      sentDown * recordBytes
    );
  }

  return NULL;
}

// Resolve local and remote candidates for owned cells
void* SugarChunk::resolveConflictsLocal(
  void* /*argument*/
)
{
  recvCandidatesFromAbove_.clear();
  recvCandidatesFromBelow_.clear();

  // Read requests from previous chunk that target this chunk
  if (previousChunk_ >= 0) {
    int offset[1] = {
        -1
    };
    void* neighbor = getOutMessage(
      1,
      offset
    );

    if (neighbor) {
      int* header = static_cast<int*>(
        neighbor
        );
      const int sentUp = header[1];
      const int sentDown = header[2];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES +
        sentUp * sizeof(
          MoveCandidateRec
          );

      recvCandidatesFromAbove_.assign(
        sentDown,
        MoveCandidateRec{}
      );

      if (sentDown > 0) {
        std::memcpy(
          recvCandidatesFromAbove_.data(),
          base,
          sentDown * sizeof(
            MoveCandidateRec
            )
        );
      }
    }
  }

  // Read requests from next chunk that target this chunk
  if (nextChunk_ >= 0) {
    int offset[1] = {
        1
    };
    void* neighbor = getOutMessage(
      1,
      offset
    );

    if (neighbor) {
      int* header = static_cast<int*>(
        neighbor
        );
      const int sentUp = header[1];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES;

      recvCandidatesFromBelow_.assign(
        sentUp,
        MoveCandidateRec{}
      );

      if (sentUp > 0) {
        std::memcpy(
          recvCandidatesFromBelow_.data(),
          base,
          sentUp * sizeof(
            MoveCandidateRec
            )
        );
      }
    }
  }

  currentStamp_++;

  if (currentStamp_ == INT_MAX) {
    std::fill(
      cellStamp_.begin(),
      cellStamp_.end(),
      0
    );
    currentStamp_ = 1;
  }

  const int localAgentCount = static_cast<int>(
    localAgents_.size()
    );

  moveWon_.assign(
    static_cast<size_t>(
      localAgentCount
      ),
    0
  );

  // Add local candidates for owned target cells
  for (int i = 0; i < localAgentCount; ++i) {
    const int targetRow = targetCells_[i].first;
    const int targetCol = targetCells_[i].second;

    if (
      !ownsRow(
        targetRow
      )
      ) {
      continue;
    }

    const int localIndex =
      (targetRow - localStartRow_) * gridWidth_ + targetCol;

    touchConflictCell(
      localIndex
    );

    if (
      localAgents_[i].id <
      bestAgentId_[
        static_cast<size_t>(
          localIndex
          )
      ]
      ) {
      bestLocalAgentIndex_[
        static_cast<size_t>(
          localIndex
          )
      ] = i;
      bestAgentId_[
        static_cast<size_t>(
          localIndex
          )
      ] = localAgents_[i].id;
    }
  }

  // Merge remote candidates using same lowest id conflict rule
  auto mergeRemote = [&](
    const std::vector<MoveCandidateRec>& candidates
    ) {
      for (const auto& record : candidates) {
        if (
          !ownsRow(
            record.targetRow
          )
          ) {
          continue;
        }

        const int localIndex =
          (record.targetRow - localStartRow_) * gridWidth_ +
          record.targetCol;

        touchConflictCell(
          localIndex
        );

        if (
          record.agentId <
          bestAgentId_[
            static_cast<size_t>(
              localIndex
              )
          ]
          ) {
          bestLocalAgentIndex_[
            static_cast<size_t>(
              localIndex
              )
          ] = -1;
          bestAgentId_[
            static_cast<size_t>(
              localIndex
              )
          ] = record.agentId;
        }
      }
    };

  mergeRemote(
    recvCandidatesFromAbove_
  );
  mergeRemote(
    recvCandidatesFromBelow_
  );

  // Mark local winners after all candidates are merged
  for (int i = 0; i < localAgentCount; ++i) {
    const int targetRow = targetCells_[i].first;
    const int targetCol = targetCells_[i].second;

    if (
      !ownsRow(
        targetRow
      )
      ) {
      continue;
    }

    const int localIndex =
      (targetRow - localStartRow_) * gridWidth_ + targetCol;

    if (
      cellStamp_[
        static_cast<size_t>(
          localIndex
          )
      ] == currentStamp_ &&
      bestLocalAgentIndex_[
        static_cast<size_t>(
          localIndex
          )
      ] == i
          ) {
      moveWon_[
        static_cast<size_t>(
          i
          )
      ] = 1;
    }
  }

  winsForFromAbove_.assign(
    recvCandidatesFromAbove_.size(),
    -1
  );
  winsForFromBelow_.assign(
    recvCandidatesFromBelow_.size(),
    -1
  );

  // Prepare winner replies for previous chunk
  for (size_t i = 0; i < recvCandidatesFromAbove_.size(); ++i) {
    const auto& record = recvCandidatesFromAbove_[i];

    if (
      !ownsRow(
        record.targetRow
      )
      ) {
      winsForFromAbove_[i] = -1;

      continue;
    }

    const int localIndex =
      (record.targetRow - localStartRow_) * gridWidth_ +
      record.targetCol;

    if (
      cellStamp_[
        static_cast<size_t>(
          localIndex
          )
      ] == currentStamp_ &&
      bestAgentId_[
        static_cast<size_t>(
          localIndex
          )
      ] == record.agentId
          ) {
      winsForFromAbove_[i] = record.agentId;
    }
  }

  // Prepare winner replies for next chunk
  for (size_t i = 0; i < recvCandidatesFromBelow_.size(); ++i) {
    const auto& record = recvCandidatesFromBelow_[i];

    if (
      !ownsRow(
        record.targetRow
      )
      ) {
      winsForFromBelow_[i] = -1;

      continue;
    }

    const int localIndex =
      (record.targetRow - localStartRow_) * gridWidth_ +
      record.targetCol;

    if (
      cellStamp_[
        static_cast<size_t>(
          localIndex
          )
      ] == currentStamp_ &&
      bestAgentId_[
        static_cast<size_t>(
          localIndex
          )
      ] == record.agentId
          ) {
      winsForFromBelow_[i] = record.agentId;
    }
  }

  return NULL;
}

// Pack conflict results for neighboring chunks
void* SugarChunk::packWinnerReplies(
  void* /*argument*/
)
{
  resetOutMessageHeader(
    PHASE_WINNER_REPLIES
  );

  int* header = static_cast<int*>(
    outMessage
    );
  const int aboveCount = static_cast<int>(
    winsForFromAbove_.size()
    );
  const int belowCount = static_cast<int>(
    winsForFromBelow_.size()
    );

  header[1] = aboveCount;
  header[2] = belowCount;

  char* base = static_cast<char*>(
    outMessage
    ) + HEADER_BYTES;

  // Above replies are stored before below replies
  if (aboveCount > 0) {
    std::memcpy(
      base,
      winsForFromAbove_.data(),
      aboveCount * sizeof(
        int
        )
    );
  }

  if (belowCount > 0) {
    std::memcpy(
      base + aboveCount * sizeof(
        int
        ),
      winsForFromBelow_.data(),
      belowCount * sizeof(
        int
        )
    );
  }

  return NULL;
}

// Apply winning moves and harvest local sugar
void* SugarChunk::applyMovesAndTakeSugar(
  void* /*argument*/
)
{
  std::vector<int> repliesForSentAbove(
    sentCandidatesAbove_.size(),
    -1
  );
  std::vector<int> repliesForSentBelow(
    sentCandidatesBelow_.size(),
    -1
  );

  // Read replies from previous chunk for candidates sent upward
  if (previousChunk_ >= 0) {
    int offset[1] = {
        -1
    };
    void* neighbor = getOutMessage(
      1,
      offset
    );

    if (neighbor) {
      int* header = static_cast<int*>(
        neighbor
        );
      const int aboveCount = header[1];
      const int belowCount = header[2];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES +
        aboveCount * sizeof(
          int
          );
      const int count = std::min(
        belowCount,
        static_cast<int>(
          repliesForSentAbove.size()
          )
      );

      if (count > 0) {
        std::memcpy(
          repliesForSentAbove.data(),
          base,
          count * sizeof(
            int
            )
        );
      }
    }
  }

  // Read replies from next chunk for candidates sent downward
  if (nextChunk_ >= 0) {
    int offset[1] = {
        1
    };
    void* neighbor = getOutMessage(
      1,
      offset
    );

    if (neighbor) {
      int* header = static_cast<int*>(
        neighbor
        );
      const int aboveCount = header[1];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES;
      const int count = std::min(
        aboveCount,
        static_cast<int>(
          repliesForSentBelow.size()
          )
      );

      if (count > 0) {
        std::memcpy(
          repliesForSentBelow.data(),
          base,
          count * sizeof(
            int
            )
        );
      }
    }
  }

  // Convert remote replies into local move wins
  for (size_t i = 0; i < sentCandidatesAbove_.size(); ++i) {
    if (
      repliesForSentAbove[i] == sentCandidatesAbove_[i].agentId &&
      i < sentCandidateIdxAbove_.size()
      ) {
      moveWon_[
        static_cast<size_t>(
          sentCandidateIdxAbove_[i]
          )
      ] = 1;
    }
  }

  for (size_t i = 0; i < sentCandidatesBelow_.size(); ++i) {
    if (
      repliesForSentBelow[i] == sentCandidatesBelow_[i].agentId &&
      i < sentCandidateIdxBelow_.size()
      ) {
      moveWon_[
        static_cast<size_t>(
          sentCandidateIdxBelow_[i]
          )
      ] = 1;
    }
  }

  const int localAgentCount = static_cast<int>(
    localAgents_.size()
    );

  // Move winning local agents and harvest local target sugar
  for (int i = 0; i < localAgentCount; ++i) {
    if (!moveWon_[i]) {
      continue;
    }

    AgentRec& agent = localAgents_[i];
    const int targetRow = targetCells_[i].first;
    const int targetCol = targetCells_[i].second;

    if (
      ownsRow(
        agent.row
      )
      ) {
      clearOccupantAt(
        agent.row,
        agent.col
      );
    }

    agent.row = targetRow;
    agent.col = targetCol;

    if (
      ownsRow(
        targetRow
      )
      ) {
      const int localIndex =
        (targetRow - localStartRow_) * gridWidth_ + targetCol;
      const int harvested = localSugar_[localIndex];

      agent.wealth += harvested;

      if (harvested > 0) {
        localSugar_[localIndex] = 0;

        markSugarGrowing(
          localIndex
        );
      }

      occupantIds_[localIndex] = agent.id;
    }
  }

  return NULL;
}

// Compute max agent records that fit in each message half
static int agentsBufferCapacity(
  int outMessageBytes
)
{
  return (
    outMessageBytes - HEADER_BYTES
    ) / static_cast<int>(
      sizeof(
        AgentRec
        )
      ) / 2;
}

// Pack agents that moved to neighboring chunks
void* SugarChunk::packMigrations(
  void* /*argument*/
)
{
  pendingMigrationsAbove_.clear();
  pendingMigrationsBelow_.clear();

  std::vector<AgentRec> keep;
  keep.reserve(
    localAgents_.size()
  );

  // Remove nonlocal agents and queue them for neighbor delivery
  for (auto& agent : localAgents_) {
    if (
      ownsRow(
        agent.row
      )
      ) {
      keep.push_back(
        agent
      );

      continue;
    }

    if (agent.row < localStartRow_) {
      pendingMigrationsAbove_.push_back(
        agent
      );
    }
    else {
      pendingMigrationsBelow_.push_back(
        agent
      );
    }
  }

  localAgents_ = std::move(
    keep
  );

  rebuildOccupantIndex();

  resetOutMessageHeader(
    PHASE_MIGRATIONS
  );

  int* header = static_cast<int*>(
    outMessage
    );
  const int aboveCount = static_cast<int>(
    pendingMigrationsAbove_.size()
    );
  const int belowCount = static_cast<int>(
    pendingMigrationsBelow_.size()
    );

  header[1] = aboveCount;
  header[2] = belowCount;

  const int capacity = agentsBufferCapacity(
    outMessageBytes_
  );
  (void)capacity;

  char* base = static_cast<char*>(
    outMessage
    ) + HEADER_BYTES;

  // Above migrations are stored before below migrations
  if (aboveCount > 0) {
    std::memcpy(
      base,
      pendingMigrationsAbove_.data(),
      aboveCount * sizeof(
        AgentRec
        )
    );
  }

  if (belowCount > 0) {
    std::memcpy(
      base + aboveCount * sizeof(
        AgentRec
        ),
      pendingMigrationsBelow_.data(),
      belowCount * sizeof(
        AgentRec
        )
    );
  }

  return NULL;
}

// Integrate agents sent from neighboring chunks
void* SugarChunk::integrateMigrations(
  void* /*argument*/
)
{
  auto readMigrations = [&](
    int directionOffset,
    int whichCount
    ) -> std::vector<AgentRec> {
      std::vector<AgentRec> result;
      int offset[1] = {
          directionOffset
      };
      void* neighbor = getOutMessage(
        1,
        offset
      );

      if (!neighbor) {
        return result;
      }

      int* header = static_cast<int*>(
        neighbor
        );
      const int aboveCount = header[1];
      const int belowCount = header[2];
      char* base = static_cast<char*>(
        neighbor
        ) + HEADER_BYTES;

      int take;
      int sourceOffset;

      // Previous chunk sends to us through its below half
      if (directionOffset < 0) {
        take = belowCount;
        sourceOffset = aboveCount * sizeof(
          AgentRec
          );
      }
      else {
        take = aboveCount;
        sourceOffset = 0;
      }

      if (take > 0) {
        result.resize(
          take
        );
        std::memcpy(
          result.data(),
          base + sourceOffset,
          take * sizeof(
            AgentRec
            )
        );
      }

      (void)whichCount;

      return result;
    };

  std::vector<AgentRec> fromAbove = readMigrations(
    -1,
    0
  );
  std::vector<AgentRec> fromBelow = readMigrations(
    1,
    0
  );

  // Add agents received from previous chunk
  for (auto& agent : fromAbove) {
    if (
      !ownsRow(
        agent.row
      )
      ) {
      continue;
    }

    const int localIndex =
      (agent.row - localStartRow_) * gridWidth_ + agent.col;
    const int harvested = localSugar_[localIndex];

    agent.wealth += harvested;

    if (harvested > 0) {
      localSugar_[localIndex] = 0;

      markSugarGrowing(
        localIndex
      );
    }

    occupantIds_[localIndex] = agent.id;

    localAgents_.push_back(
      agent
    );
  }

  // Add agents received from next chunk
  for (auto& agent : fromBelow) {
    if (
      !ownsRow(
        agent.row
      )
      ) {
      continue;
    }

    const int localIndex =
      (agent.row - localStartRow_) * gridWidth_ + agent.col;
    const int harvested = localSugar_[localIndex];

    agent.wealth += harvested;

    if (harvested > 0) {
      localSugar_[localIndex] = 0;

      markSugarGrowing(
        localIndex
      );
    }

    occupantIds_[localIndex] = agent.id;

    localAgents_.push_back(
      agent
    );
  }

  return NULL;
}

// Harvest sugar from each agent's current cell
void* SugarChunk::takeSugarAtCurrent(
  void* /*argument*/
)
{
  for (auto& agent : localAgents_) {
    if (
      !ownsRow(
        agent.row
      )
      ) {
      continue;
    }

    const int localIndex =
      (agent.row - localStartRow_) * gridWidth_ + agent.col;

    if (occupantIds_[localIndex] != agent.id) {
      continue;
    }

    const int harvested = localSugar_[localIndex];

    if (harvested > 0) {
      agent.wealth += harvested;
      localSugar_[localIndex] = 0;

      markSugarGrowing(
        localIndex
      );
    }
  }

  return NULL;
}

// Apply metabolism and remove dead agents
void* SugarChunk::metabolismAndDeath(
  void* /*argument*/
)
{
  std::vector<AgentRec> keep;
  keep.reserve(
    localAgents_.size()
  );

  for (auto& agent : localAgents_) {
    agent.wealth -= agent.metabolism;

    const bool alive = agent.wealth > 0;

    if (alive) {
      keep.push_back(
        agent
      );
    }
    else {
      if (
        ownsRow(
          agent.row
        )
        ) {
        clearOccupantAt(
          agent.row,
          agent.col
        );
      }
    }
  }

  localAgents_ = std::move(
    keep
  );

  return NULL;
}

// Keep agent order stable for reports
void* SugarChunk::sortLocalAgents(
  void* /*argument*/
)
{
  std::sort(
    localAgents_.begin(),
    localAgents_.end(),
    [](
      const AgentRec& left,
      const AgentRec& right
      ) {
        return left.id < right.id;
    }
  );

  return NULL;
}

// Set timestep from driver
void* SugarChunk::setTimestep(
  void* argument
)
{
  if (!argument) {
    return NULL;
  }

  SugarChunkStepArg* stepArg = static_cast<SugarChunkStepArg*>(
    argument
    );

  currentStep_ = stepArg->step;

  return NULL;
}

// Advance timestep for compound MASS calls
void* SugarChunk::advanceTimestep(
  void* /*argument*/
)
{
  if (!compoundStepPrimed_) {
    compoundStepPrimed_ = true;
    currentStep_ = 0;
  }
  else {
    currentStep_++;
  }

  return NULL;
}

// Report local totals for final reduction
void* SugarChunk::reportStats(
  void* /*argument*/
)
{
  reportStatsBuf_.localAgentCount = static_cast<int>(
    localAgents_.size()
    );
  reportStatsBuf_.localWealth = 0;

  for (const auto& agent : localAgents_) {
    reportStatsBuf_.localWealth += agent.wealth;
  }

  reportStatsBuf_.localSugar = 0;

  for (int value : localSugar_) {
    reportStatsBuf_.localSugar += value;
  }

  reportStatsBuf_.localMaxAgentId = -1;

  for (const auto& agent : localAgents_) {
    if (agent.id > reportStatsBuf_.localMaxAgentId) {
      reportStatsBuf_.localMaxAgentId = agent.id;
    }
  }

  return &reportStatsBuf_;
}

// Report local agent count before agent data transfer
void* SugarChunk::reportAgentsCount(
  void* /*argument*/
)
{
  reportCountBuf_ = static_cast<int>(
    localAgents_.size()
    );

  return &reportCountBuf_;
}

// Report local agents into fixed driver buffer
void* SugarChunk::reportAgents(
  void* /*argument*/
)
{
  if (
    static_cast<int>(
      reportAgentsBuf_.size()
      ) < reportAgentsCapacityBytes_
    ) {
    reportAgentsBuf_.assign(
      reportAgentsCapacityBytes_,
      0
    );
  }

  const int agentCount = static_cast<int>(
    localAgents_.size()
    );
  const int capacityAgents = (
    reportAgentsCapacityBytes_ - static_cast<int>(
      sizeof(
        int
        )
      )
    ) / static_cast<int>(
      sizeof(
        AgentRec
        )
      );
  char* buffer = reportAgentsBuf_.data();

  if (agentCount > capacityAgents) {
    std::fprintf(
      stderr,
      "ERROR: SugarChunk %d reportAgents_ overflow: have %d agents but "
      "buffer holds only %d. The deterministic V2 broadcast cannot "
      "silently truncate; the driver will abort.\n",
      chunkIndex_,
      agentCount,
      capacityAgents
    );

    *reinterpret_cast<int*>(
      buffer
      ) = -1;

    return buffer;
  }

  // First int stores count, followed by raw AgentRec data
  *reinterpret_cast<int*>(
    buffer
    ) = agentCount;

  if (agentCount > 0) {
    std::memcpy(
      buffer + sizeof(
        int
        ),
      localAgents_.data(),
      agentCount * sizeof(
        AgentRec
        )
    );
  }

  return buffer;
}

// Report local sugar rows into fixed driver buffer
void* SugarChunk::reportSugar(
  void* /*argument*/
)
{
  if (
    static_cast<int>(
      reportSugarBuf_.size()
      ) < reportSugarCapacityInts_
    ) {
    reportSugarBuf_.assign(
      reportSugarCapacityInts_,
      0
    );
  }
  else {
    std::fill(
      reportSugarBuf_.begin(),
      reportSugarBuf_.end(),
      0
    );
  }

  // First 2 ints identify the owned row block
  reportSugarBuf_[0] = localStartRow_;
  reportSugarBuf_[1] = localRowCount_;

  int total = localRowCount_ * gridWidth_;
  const int maxPayload = reportSugarCapacityInts_ - 2;

  if (total > maxPayload) {
    total = maxPayload;
  }

  if (total > 0) {
    std::memcpy(
      reportSugarBuf_.data() + 2,
      localSugar_.data(),
      total * sizeof(
        int
        )
    );
  }

  return reportSugarBuf_.data();
}