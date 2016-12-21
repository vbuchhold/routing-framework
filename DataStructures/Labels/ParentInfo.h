#pragma once

// This enum specifies the different kinds of parent information that can be collected by Dijkstra.
enum class ParentInfo {
  NO_PARENT_INFO,
  PARENT_VERTICES_ONLY,
  FULL_PARENT_INFO,
};
