package rw

import (
	"bytes"
	"encoding/binary"
	"math"
)

type ReaderWriter struct {
	order   binary.ByteOrder
	dataBuf []byte
	rw      bytes.Buffer
}

func NewNavMeshDataBinWriter() *ReaderWriter {
	return &ReaderWriter{order: binary.LittleEndian, dataBuf: make([]byte, 8)}
}

func NewNavMeshDataBinReader(data []byte) *ReaderWriter {
	d := &ReaderWriter{order: binary.LittleEndian, dataBuf: make([]byte, 8)}
	d.rw.Write(data)
	return d
}

func (w *ReaderWriter) ReadInt16() int16 {
	return int16(w.ReadUInt16())
}

func (w *ReaderWriter) ReadInt16s(value []int16) {
	for i := range value {
		value[i] = w.ReadInt16()
	}
}

func (w *ReaderWriter) ReadUInt16s(value []uint16) {
	for i := range value {
		value[i] = w.ReadUInt16()
	}
}

func (w *ReaderWriter) ReadUInt16() uint16 {
	_, err := w.rw.Read(w.dataBuf[:2])
	if err != nil {
		panic(err)
	}
	return w.order.Uint16(w.dataBuf[:2])
}

func (w *ReaderWriter) ReadInt8s(value []int8) {
	for i := range value {
		value[i] = w.ReadInt8()
	}
}

func (w *ReaderWriter) ReadUInt8s(value []uint8) {
	for i := range value {
		value[i] = w.ReadUInt8()
	}
}

func (w *ReaderWriter) ReadInt8() int8 {
	return int8(w.ReadUInt8())
}

func (w *ReaderWriter) ReadUInt8() uint8 {
	res, err := w.rw.ReadByte()
	if err != nil {
		panic(err)
	}
	return res
}

func (w *ReaderWriter) ReadInt32s(value []int32) {
	for i := range value {
		value[i] = w.ReadInt32()
	}
}

func (w *ReaderWriter) ReadUInt32s(value []uint32) {
	for i := range value {
		value[i] = w.ReadUInt32()
	}
}
func (w *ReaderWriter) ReadInt32() int32 {
	return int32(w.ReadUInt32())
}

func (w *ReaderWriter) ReadUInt32() uint32 {
	_, err := w.rw.Read(w.dataBuf[:4])
	if err != nil {
		panic(err)
	}
	return w.order.Uint32(w.dataBuf[:4])
}
func (w *ReaderWriter) ReadFloat32s(value []float32) {
	for i := range value {
		value[i] = w.ReadFloat32()
	}
}
func (w *ReaderWriter) ReadFloat32() float32 {
	return math.Float32frombits(w.ReadUInt32())
}

func (w *ReaderWriter) WriteInt16(v interface{}) {
	switch value := v.(type) {
	case int16:
		w.order.PutUint16(w.dataBuf, uint16(value))
	case uint16:
		w.order.PutUint16(w.dataBuf, value)
	default:
		panic("not impl")
	}
	w.rw.Write(w.dataBuf[:2])
}

func (w *ReaderWriter) WriteInt16s(v interface{}) {
	switch value := v.(type) {
	case []int16:
		for _, tmp := range value {
			w.WriteInt16(tmp)
		}
	case []uint16:
		for _, tmp := range value {
			w.WriteInt16(tmp)
		}
	default:
		panic("not impl")
	}
	w.rw.Write(w.dataBuf[:2])
}

func (w *ReaderWriter) WriteInt8(v interface{}) {
	switch value := v.(type) {
	case int8:
		w.rw.WriteByte(byte(value))
	case uint8:
		w.rw.WriteByte(value)
	default:
		panic("not impl")
	}
}

func (w *ReaderWriter) WriteInt8s(v interface{}) {
	switch value := v.(type) {
	case []int8:
		for _, tmp := range value {
			w.WriteInt8(tmp)
		}
	case []uint8:
		for _, tmp := range value {
			w.WriteInt8(tmp)
		}
	default:
		panic("not impl")
	}
}

func (w *ReaderWriter) WriteInt32(v interface{}) {
	switch value := v.(type) {
	case int32:
		w.order.PutUint32(w.dataBuf, uint32(value))
	case int:
		w.order.PutUint32(w.dataBuf, uint32(value))
	case uint32:
		w.order.PutUint32(w.dataBuf, value)
	default:
		panic("not impl")
	}
	w.rw.Write(w.dataBuf[:4])
}

func (w *ReaderWriter) WriteInt32s(v interface{}) {
	switch value := v.(type) {
	case []int32:
		for _, tmp := range value {
			w.WriteInt32(tmp)
		}
	case []int:
		for _, tmp := range value {
			w.WriteInt32(tmp)
		}
	case []uint32:
		for _, tmp := range value {
			w.WriteInt32(tmp)
		}
	default:
		panic("not impl")
	}

}

func (w *ReaderWriter) WriteFloat32(v interface{}) {
	switch value := v.(type) {
	case float32:
		w.order.PutUint32(w.dataBuf, math.Float32bits(value))
	case float64:
		w.order.PutUint32(w.dataBuf, math.Float32bits(float32(value)))
	default:
		panic("not impl")
	}
	w.rw.Write(w.dataBuf[:4])
}

func (w *ReaderWriter) WriteFloat32s(v interface{}) {
	switch value := v.(type) {
	case []float32:
		for _, tmp := range value {
			w.WriteFloat32(tmp)
		}
	case []float64:
		for _, tmp := range value {
			w.WriteFloat32(float32(tmp))
		}
	default:
		panic("not impl")
	}
}
func (w *ReaderWriter) Skip(size int) {
	w.rw.Next(size)
}
func (w *ReaderWriter) GetWriteBytes() (res []byte) {
	res = w.rw.Bytes()
	return res
}

func (w *ReaderWriter) PadZero(n int) {
	for i := 0; i < n; i++ {
		w.rw.WriteByte(0)
	}
}
func (w *ReaderWriter) ChangeOrder(order binary.ByteOrder) {
	w.order = order
}

func (w *ReaderWriter) Size() int {
	return w.rw.Len()
}
