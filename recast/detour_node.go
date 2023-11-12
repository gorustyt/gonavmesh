package recast

import (
	"container/heap"
)

const (
	DT_NODE_OPEN            = 0x01
	DT_NODE_CLOSED          = 0x02
	DT_NODE_PARENT_DETACHED = 0x04 // parent of the node is not adjacent. Found using raycast.
)

type dtNodeIndex int

const (
	DT_NODE_PARENT_BITS    = 24
	DT_NODE_STATE_BITS     = 2
	DT_MAX_STATES_PER_NODE = 1 << DT_NODE_STATE_BITS // number of extra states per node. See dtNode::state
)

type dtNode struct {
	pos    [3]float64 ///< Position of the node.
	cost   float64    ///< Cost from previous node to current node.
	total  float64    ///< Cost up to the node.
	pidx   *dtNode    ///<  parent node.
	state  int        ///< extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE
	flags  int        ///< Node flags. A combination of dtNodeFlags.
	id     dtPolyRef  ///< Polygon ref the node corresponds to.
	_index int        //堆中移除和更新对象用
}

func (node *dtNode) SetIndex(index int) {
	node._index = index
}

func (node *dtNode) GetIndex() int {
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
func (q *nodeQueue[T]) Poll() T { return heap.Pop(q) } //从堆顶弹出一个元素
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
	q.data = append(q.data, x)
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
