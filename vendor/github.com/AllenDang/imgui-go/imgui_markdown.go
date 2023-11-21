package imgui

// #include "imgui_markdown_wrapper.h"
import "C"

type MarkdownHeaderData struct {
	Font         Font
	HasSeparator bool
}

type MarkdownImageData struct {
	TextureID              *TextureID
	Scale                  bool
	Size                   Vec2
	UseLinkCallback        bool
	Uv0, Uv1               Vec2
	TintColor, BorderColor Vec4
}

// markdownImageCallbackCache stores user-definied image loader
// it is only way to share it with C callback
var markdownImageCallbackCache func(url string) MarkdownImageData

// markdownImageCache stores markdown image data
// TODO: meybe it should be done on client-side?...
var markdownImageCache map[string]*MarkdownImageData

func init() {
	markdownImageCache = make(map[string]*MarkdownImageData)
}

// Markdown implements imgui_markdown.h
// NOTE about arguments:
// - data is pointer to markdown text data
// - linkCB is callback called when link is clicked (it should most likely open link in a web browser)
// - imageCB is expected to load an image at `path` and return POINTER to its texture with some other
//   stats. BE AWARE that imageCB MUST NOT pause goroutine. It could e.g.:
//   - create variable with TextureID = 0
//   - invoge a new goroutine and load the texture there and set pointers value
//   - return pointer to declared variable
// - headers are headers formatting data. Note, that first index of slice will be applied
//   to top-level (H1), second for H2 and so on.
func Markdown(data *string, linkCB func(s string), imageCB func(path string) MarkdownImageData, headers []MarkdownHeaderData) {
	// share imageCB with C callback (goMarkdownImageCallback)
	markdownImageCallbackCache = imageCB

	state := newInputTextState(*data, nil)
	defer func() {
		*data = state.buf.toGo()
		state.release()
	}()

	// prepare headers for C
	cHeaders := []C.iggMarkdownHeaderData{}
	if headers != nil {
		for _, data := range headers {
			cHeaders = append(cHeaders,
				C.iggMarkdownHeaderData{
					font:      data.Font.handle(),
					separator: castBool(data.HasSeparator),
				},
			)
		}
	}

	var cHeadersPtr *C.iggMarkdownHeaderData
	if len(cHeaders) > 0 {
		// this trick allows to pass go slice into C
		// documentation: https://coderwall.com/p/m_ma7q/pass-go-slices-as-c-array-parameters
		cHeadersPtr = &cHeaders[0]
	}

	linkData := C.iggMarkdown(
		(*C.char)(state.buf.ptr),
		cHeadersPtr, (C.int)(len(cHeaders)),
	)

	// Read link callback
	s := C.GoString(linkData.link)
	s = s[:int(linkData.link_len)]
	if s != "" {
		linkCB(s)
	}
}

// goMarkdownImageCallback is exported to C callback for loading markdown images.
// in short, it calls user-definied cached in markdownImageCallbackCache function.
//export goMarkdownImageCallback
func goMarkdownImageCallback(data C.iggMarkdownLinkCallbackData) (result C.iggMarkdownImageData) {
	if markdownImageCallbackCache == nil {
		return result
	}

	path := C.GoString(data.link)
	path = path[:int(data.link_len)]

	// it calls user-definied function only at first time when this is called.
	if _, found := markdownImageCache[path]; !found {
		markdownImageCache[path] = &MarkdownImageData{}
		go func() {
			d := markdownImageCallbackCache(path)
			*markdownImageCache[path] = d
		}()
	}

	d := markdownImageCache[path]

	// check if texture id isn't nil
	if d.TextureID == nil {
		return result
	}

	sizeArg, _ := d.Size.wrapped()
	uv0, _ := d.Uv0.wrapped()
	uv1, _ := d.Uv1.wrapped()
	tintColor, _ := d.TintColor.wrapped()
	borderColor, _ := d.BorderColor.wrapped()

	result = C.iggMarkdownImageData{
		texture:         d.TextureID.handle(),
		shouldScale:     castBool(d.Scale),
		size:            *sizeArg,
		useLinkCallback: castBool(d.UseLinkCallback),
		uv0:             *uv0,
		uv1:             *uv1,
		tintColor:       *tintColor,
		borderColor:     *borderColor,
	}

	// return to C
	return result
}
