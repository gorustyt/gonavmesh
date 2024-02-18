package detour

import (
	"container/heap"
	"github.com/gorustyt/gonavmesh/common"
)

const (
	DT_NODE_OPEN            = 0x01
	DT_NODE_CLOSED          = 0x02
	DT_NODE_PARENT_DETACHED = 0x04 // parent of the node is not adjacent. Found using raycast.
)

type DtNodeIndex int

const (
	DT_NODE_PARENT_BITS    = 24
	DT_NODE_STATE_BITS     = 2
	DT_MAX_STATES_PER_NODE = 1 << DT_NODE_STATE_BITS // number of extra states per node. See DtNode::state
)

type DtNode struct {
	Pos    []float32 ///< Position of the node.
	Cost   float32   ///< Cost from previous node to current node.
	Total  float32   ///< Cost up to the node.
	Pidx   uint32    ///<  parent node.
	State  uint32    ///< extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE
	Flags  uint32    ///< Node flags. A combination of DtNodeFlags.
	Id     DtPolyRef ///< Polygon ref the node corresponds to.
	_index int       //堆中移除和更新对象用
}

func NewDtNode() *DtNode {
	return &DtNode{
		Pidx:  DT_NODE_PARENT_BITS,
		State: DT_NODE_STATE_BITS,
		Flags: 3,
		Pos:   make([]float32, 3),
	}
}

func (node *DtNode) SetIndex(index int) {
	node._index = index
}

func (node *DtNode) GetIndex() int {
	return node._index
}

type NodeQueueIndex interface {
	SetIndex(index int)
	GetIndex() int
}
type NodeQueue[T any] interface {
	Peek() T         //查看堆顶，不会移除元素
	Poll() T         //从堆顶弹出一个元素
	Update(any) bool //更新元素
	Remove(any) bool //移除一个元素
	Offer(T)         //插入一个元素
	Reset()
	Empty() bool
}

// 优先级队列
type nodeQueue[T any] struct {
	data []T
	less func(t1, t2 T) bool
}

func NewNodeQueue[T any](less func(t1, t2 T) bool) NodeQueue[T] {
	q := &nodeQueue[T]{less: less}
	heap.Init(q)
	return q
}
func (q *nodeQueue[T]) Reset() {
	q.data = []T{}
}

// 查看堆顶
func (q *nodeQueue[T]) Peek() T {
	return q.data[0]
}
func (q *nodeQueue[T]) Poll() T { return heap.Pop(q).(T) } //从堆顶弹出一个元素
// 更新元素
func (q *nodeQueue[T]) Update(value any) bool {

	if v, ok := value.(NodeQueueIndex); ok {
		heap.Fix(q, v.GetIndex())
		return true
	}
	return false

}

func (q *nodeQueue[T]) Remove(value any) bool {
	if v, ok := value.(NodeQueueIndex); ok {
		heap.Remove(q, v.GetIndex())
		return true
	}
	return false

}                                     //移除一个元素
func (q *nodeQueue[T]) Offer(value T) { heap.Push(q, value) } //插入一个元素

func (q *nodeQueue[T]) Push(x any) {
	q.data = append(q.data, x.(T))
	if v, ok := x.(NodeQueueIndex); ok {
		v.SetIndex(len(q.data) - 1)
	}
}
func (q *nodeQueue[T]) Pop() (res any) {
	res = q.data[0]
	if v, ok := res.(NodeQueueIndex); ok {
		v.SetIndex(0)
	}
	return res
}
func (q *nodeQueue[T]) Len() int {
	return len(q.data)
}
func (q *nodeQueue[T]) Empty() bool {
	return q.Len() == 0
}
func (q *nodeQueue[T]) Less(i, j int) bool { return q.less(q.data[i], q.data[j]) }
func (q *nodeQueue[T]) Swap(i, j int) {
	var vi any = q.data[i]
	var vj any = q.data[j]
	if v, ok := vi.(NodeQueueIndex); ok {
		v.SetIndex(j)
	}
	if v, ok := vj.(NodeQueueIndex); ok {
		v.SetIndex(i)
	}
	q.data[i], q.data[j] = q.data[j], q.data[i]

}

func dtHashRef(a DtPolyRef) int32 {
	a += ^(a << 15)
	a ^= (a >> 10)
	a += (a << 3)
	a ^= (a >> 6)
	a += ^(a << 11)
	a ^= (a >> 16)
	return int32(a)
}

const (
	DT_NULL_IDX = DtNodeIndex(^0)
)

type DtNodePool struct {
	m_nodes     []*DtNode
	m_first     []DtNodeIndex
	m_next      []DtNodeIndex
	m_maxNodes  int32
	m_hashSize  int32
	m_nodeCount int32
}

func NewDtNodePool(maxNodes, hashSize int32) *DtNodePool {
	p := &DtNodePool{
		m_maxNodes: maxNodes,
		m_hashSize: hashSize,
	}
	common.AssertTrue(common.NextPow2(uint32(p.m_hashSize)) == uint32(p.m_hashSize))
	// pidx is special as 0 means "none" and 1 is the first node. For that reason
	// we have 1 fewer nodes available than the number of values it can contain.
	common.AssertTrue(p.m_maxNodes > 0 && int(p.m_maxNodes) <= int(DT_NULL_IDX) && p.m_maxNodes <= (1<<DT_NODE_PARENT_BITS)-1)
	p.m_nodes = make([]*DtNode, p.m_maxNodes)
	p.m_next = make([]DtNodeIndex, p.m_maxNodes)
	p.m_first = make([]DtNodeIndex, hashSize)
	for i := range p.m_next {
		p.m_next[i] = 0xff
	}
	for i := range p.m_first {
		p.m_first[i] = 0xff
	}
	return p
}
func (p *DtNodePool) GetNodeIdx(node *DtNode) int {
	if node == nil {
		return 0
	}
	return p.getIndexByNode(node)
}

func (p *DtNodePool) getIndexByNode(node *DtNode) int {
	for i, v := range p.m_nodes {
		if v == node {
			return i
		}
	}
	return -1
}

func (p *DtNodePool) GetNodeAtIdx(idx int) *DtNode {
	if idx == 0 {
		return nil
	}
	return p.m_nodes[idx-1]
}

func (p *DtNodePool) GetMaxNodes() int32 { return p.m_maxNodes }

func (p *DtNodePool) GetHashSize() int32              { return p.m_hashSize }
func (p *DtNodePool) GetFirst(bucket int) DtNodeIndex { return p.m_first[bucket] }
func (p *DtNodePool) GetNext(i int) DtNodeIndex       { return p.m_next[i] }
func (p *DtNodePool) GetNodeCount() int32             { return p.m_nodeCount }
func (p *DtNodePool) Clear() {
	p.m_next = make([]DtNodeIndex, len(p.m_next))
	p.m_nodeCount = 0
}

func (p *DtNodePool) FindNodes(id DtPolyRef, maxNodes int) (nodes []*DtNode, n int) {
	bucket := dtHashRef(id) & (p.m_hashSize - 1)
	i := p.m_first[bucket]
	for i != DT_NULL_IDX {
		if p.m_nodes[i].Id == id {
			if n >= maxNodes {
				return nodes, n
			}
			nodes = append(nodes, p.m_nodes[i])
			n++
		}
		i = p.m_next[i]
	}

	return nodes, n
}

func (p *DtNodePool) GetNode(id DtPolyRef, states ...uint32) *DtNode {
	state := uint32(0)
	if len(states) > 0 {
		state = states[0]
	}
	bucket := dtHashRef(id) & (p.m_hashSize - 1)
	i := p.m_first[bucket]
	for i != DT_NULL_IDX {
		if p.m_nodes[i].Id == id && p.m_nodes[i].State == state {
			return p.m_nodes[i]
		}

		i = p.m_next[i]
	}

	if p.m_nodeCount >= p.m_maxNodes {
		return nil
	}

	i = DtNodeIndex(p.m_nodeCount)
	p.m_nodeCount++

	// Init node
	node := p.m_nodes[i]
	if node == nil {
		node = NewDtNode()
		p.m_nodes[i] = node
	}
	node.Pidx = 0
	node.Cost = 0
	node.Total = 0
	node.Id = id
	node.State = state
	node.Flags = 0

	p.m_next[i] = p.m_first[bucket]
	p.m_first[bucket] = i

	return node
}

func (p *DtNodePool) FindNode1(id DtPolyRef, state uint32) *DtNode {
	bucket := dtHashRef(id) & (p.m_hashSize - 1)
	i := p.m_first[bucket]
	for i != DT_NULL_IDX {
		if p.m_nodes[i].Id == id && p.m_nodes[i].State == state {
			return p.m_nodes[i]
		}

		i = p.m_next[i]
	}
	return nil
}
