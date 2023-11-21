#include "imnodes.h"
#include "imnodesWrapper.h"
#include "imguiWrappedHeader.h"
#include "WrapperConverter.h"

void iggImNodesCreateContext()
{
  ImNodes::CreateContext();
}

void iggImNodesDestroyContext()
{
  ImNodes::DestroyContext();
}

void iggImNodesBeginNodeEditor()
{
  ImNodes::BeginNodeEditor();
}

void iggImNodesEndNodeEditor()
{
  ImNodes::EndNodeEditor();
}

void iggImNodesBeginNode(int id)
{
  ImNodes::BeginNode(id);
}

void iggImNodesEndNode()
{
  ImNodes::EndNode();
}

void iggImNodesBeginNodeTitleBar()
{
  ImNodes::BeginNodeTitleBar();
}

void iggImNodesEndNodeTitleBar()
{
  ImNodes::EndNodeTitleBar();
}

void iggImNodesBeginInputAttribute(int id)
{
  ImNodes::BeginInputAttribute(id);
}

void iggImNodesEndInputAttribute()
{
  ImNodes::EndInputAttribute();
}

void iggImNodesBeginOutputAttribute(int id)
{
  ImNodes::BeginOutputAttribute(id);
}

void iggImNodesEndOutputAttribute()
{
  ImNodes::EndOutputAttribute();
}

void iggImNodesLink(int id, int start_attribute_id, int end_attribute_id)
{
  ImNodes::Link(id, start_attribute_id, end_attribute_id);
}

IggBool iggImNodesIsLinkCreated(
    int* started_at_node_id,
    int* started_at_attribute_id,
    int* ended_at_node_id,
    int* ended_at_attribute_id,
    IggBool* created_from_snap)
{
  BoolWrapper boolArg(created_from_snap);
  return ImNodes::IsLinkCreated(started_at_node_id, started_at_attribute_id, ended_at_node_id, ended_at_attribute_id, boolArg) ? 1 : 0;
}

IggBool iggImNodesIsLinkDestroyed(int* link_id)
{
  return ImNodes::IsLinkDestroyed(link_id) ? 1 : 0;
}

void iggImNodesPushAttributeFlag(int flag)
{
  ImNodes::PushAttributeFlag(static_cast<ImNodesAttributeFlags>(flag));
}

void iggImNodesPopAttributeFlag()
{
  ImNodes::PopAttributeFlag();
}

void iggImNodesEnableDetachWithCtrlClick()
{
  ImNodesIO& io = ImNodes::GetIO();
  io.LinkDetachWithModifierClick.Modifier = &ImGui::GetIO().KeyCtrl;
}

void iggImNodesSetNodeScreenSpacePos(int node_id, const IggVec2 *screen_space_pos)
{
  Vec2Wrapper posArg(screen_space_pos);
  ImNodes::SetNodeScreenSpacePos(node_id, *posArg);
}

void iggImNodesSetNodeEditorSpacePos(int node_id, const IggVec2 *editor_space_pos)
{
  Vec2Wrapper posArg(editor_space_pos);
  ImNodes::SetNodeEditorSpacePos(node_id, *posArg);
}

void iggImNodesSetNodeGridSpacePos(int node_id, const IggVec2 *grid_pos)
{
  Vec2Wrapper posArg(grid_pos);
  ImNodes::SetNodeGridSpacePos(node_id, *posArg);
}

void iggImNodesGetNodeScreenSpacePos(const int node_id, IggVec2 *pos)
{
  exportValue(*pos, ImNodes::GetNodeScreenSpacePos(node_id));
}

void iggImNodesGetNodeEditorSpacePos(const int node_id, IggVec2 *pos)
{
  exportValue(*pos, ImNodes::GetNodeEditorSpacePos(node_id));
}

void iggImNodesGetNodeGridSpacePos(const int node_id, IggVec2 *pos)
{
  exportValue(*pos, ImNodes::GetNodeGridSpacePos(node_id));
}
