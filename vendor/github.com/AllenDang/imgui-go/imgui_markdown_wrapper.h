#pragma once

#include "imguiWrapperTypes.h"

#ifdef __cplusplus
extern "C"
{
#endif
	typedef struct {
		char* link;
		int link_len;
	} iggMarkdownLinkCallbackData;

        typedef struct {
                IggFont font;
                IggBool separator;
        } iggMarkdownHeaderData;

        typedef struct {
                IggTextureID texture;
                IggBool useLinkCallback;
                IggVec2 size;
                IggBool shouldScale; // if true, image will not be larger than available region

                IggVec2 uv0;
                IggVec2 uv1;
                IggVec4 tintColor;
                IggVec4 borderColor;
        } iggMarkdownImageData;

        extern iggMarkdownLinkCallbackData iggMarkdown(
                        char *markdown_,
                        iggMarkdownHeaderData fonts[], int numHeaderLevels
                );

        extern iggMarkdownImageData goMarkdownImageCallback(iggMarkdownLinkCallbackData);
#ifdef __cplusplus
}
#endif
