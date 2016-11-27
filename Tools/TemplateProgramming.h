#pragma once

namespace impl {

using PackExpansionT = int[];

}

// Let pattern be a statement naming one or more parameter packs. This macro executes the pattern
// for each element in the pack.
#define RUN_FORALL(pattern) impl::PackExpansionT{0, ((pattern), void(), 0)...}
