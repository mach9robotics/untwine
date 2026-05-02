#pragma once

#include <stdexcept>

namespace untwine
{

class  FatalError : public std::runtime_error
{
public:
    inline FatalError(std::string const& msg) : std::runtime_error(msg)
        {}
};

// Thrown when a just-written COPC chunk fails to round-trip through
// lazperf's decompressor at its declared point count. Signals a rare
// non-deterministic lazperf bug; callers should retry untwine.
//
// main() translates this to a distinct non-zero exit code so external
// orchestrators (e.g. the ingest pipeline) can detect chunk corruption
// specifically and retry rather than abort the run entirely.
class ChunkIntegrityError : public FatalError
{
public:
    inline ChunkIntegrityError(std::string const& msg) : FatalError(msg)
        {}
};

// Process exit code indicating that untwine detected a corrupt COPC chunk
// (declared point_count exceeds what the compressed LAZ bytes decode to).
// Chosen to be distinct from generic error exit (-1/255) so orchestrators
// can distinguish "retryable corruption" from "hard failure".
constexpr int ChunkIntegrityExitCode = 43;

} // namespace untwine
