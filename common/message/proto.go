package message

import (
	"gonavamesh/common"
	"google.golang.org/protobuf/proto"
)

func Encode(msg proto.Message) (data []byte) {
	data, err := proto.Marshal(msg)
	common.AssertTrue(err == nil)
	return data
}

func Decode(data []byte, msg proto.Message) {
	err := proto.Unmarshal(data, msg)
	common.AssertTrue(err == nil)
}
