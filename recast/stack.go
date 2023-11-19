package recast

type Stack[T any] interface {
	Pop() T
	Push(value T)
	Len() int
	Empty() bool
	Clear()
	Index(index int) T
	SetByIndex(index int, value T)
	Resize(size int, value ...T)
	Reserve(count int)
	Slice(begin, end int) []T
	Data() []T
}

func NewStack[T any](construct func() T) Stack[T] {
	return NewStackArray(construct, 0)
}

func NewStackArray[T any](construct func() T, count int) Stack[T] {
	s := &stack[T]{construct: construct}
	s.Resize(count)
	return s
}

type stack[T any] struct {
	data      []T
	construct func() T //T的构造函数
	m_cap     int
}

func (s *stack[T]) Data() []T {
	return s.data
}
func (s *stack[T]) Clear() {
	s.data = []T{}
}
func (s *stack[T]) Pop() T {
	e := s.data[s.Len()-1]
	s.data = s.data[:s.Len()-1]
	return e
}

func (s *stack[T]) Push(value T) {
	s.data = append(s.data, value)
}

func (s *stack[T]) Len() int {
	return len(s.data)
}

func (s *stack[T]) Empty() bool {
	return s.Len() == 0
}

func (s *stack[T]) Index(index int) T {
	return s.data[index]
}

// [begin.end)
func (s *stack[T]) Slice(begin, end int) []T {
	if end > s.Len() {
		end = s.Len()
	}
	return s.data[begin:end]
}
func (s *stack[T]) SetByIndex(index int, value T) {
	s.data[index] = value
}
func (s *stack[T]) Reserve(count int) {
	if count >= s.Len() {
		return
	}
	s.data = s.data[:count]
}
func (s *stack[T]) Resize(size int, value ...T) {
	m_size := s.Len()
	if size < s.Len() {
		s.data = s.data[:size]
	} else if size > m_size {
		for i := 0; i < size-m_size; i++ {
			if len(value) > 0 {
				s.Push(value[0])
			} else {
				s.Push(s.construct())
			}
		}
		if size > s.m_cap {
			s.m_cap = size
		}
	}
}
