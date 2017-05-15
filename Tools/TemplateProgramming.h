#pragma once

namespace templates {

using PackExpansion = int[];

}

// Let pattern be an expression naming one or more parameter packs. This macro executes the pattern
// for each element in the pack.
#define RUN_FORALL(pattern) templates::PackExpansion{0, ((pattern), void(), 0)...}

// Let pattern be an expression naming one or more parameter packs. This macro executes the pattern
// for each element in the pack for which cond holds.
#define RUN_IF(cond, pattern) RUN_FORALL(((cond) ? (pattern), void() : void()))

// An empty class, intended for implementing conditional inheritance.
class EmptyClass {};
