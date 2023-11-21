/*
 * =====================================================================================
 *
 *       Filename:  imgui_markdown_wrapper.cpp
 *
 *    Description:  wrapper for imgui_markdown.h
 *
 * =====================================================================================
 */

#include <WrapperConverter.h>
#include <imgui.h>
#include <imgui_markdown.h>
#include <imgui_markdown_wrapper.h>

iggMarkdownLinkCallbackData iggMarkdownLink;
static void markdownLinkCallback(ImGui::MarkdownLinkCallbackData);
static ImGui::MarkdownImageData markdownImageCallback(ImGui::MarkdownLinkCallbackData);

iggMarkdownLinkCallbackData iggMarkdown(char *markdown_, iggMarkdownHeaderData headers_[], int numHeaderLevels) {
  // clean up link cache.
  iggMarkdownLink.link = NULL;
  iggMarkdownLink.link_len = 0;

  // create imgui markdown config:
  // TODO: implement all these methods
  ImGui::MarkdownConfig mdConfig;

  mdConfig.linkCallback = markdownLinkCallback;
  for (int i = 0; i < numHeaderLevels; i++) {
    ImFont *font = reinterpret_cast<ImFont *>(headers_[i].font);
    mdConfig.headingFormats[i] = {font, (bool)(headers_[i].separator)};
  }
  // mdConfig.tooltipCallback =      NULL;
  mdConfig.imageCallback = markdownImageCallback;
  // mdConfig.linkIcon =             ICON_FA_LINK;
  // mdConfig.userData =             NULL;
  // mdConfig.formatCallback =       ExampleMarkdownFormatCallback;

  // run ImGui Markdown
  ImGui::Markdown(markdown_, strlen(markdown_), mdConfig);

  return iggMarkdownLink;
}

void markdownLinkCallback(ImGui::MarkdownLinkCallbackData data_) {
  iggMarkdownLink.link = (char *)(data_.link);
  iggMarkdownLink.link_len = data_.linkLength;
}

ImGui::MarkdownImageData markdownImageCallback(ImGui::MarkdownLinkCallbackData data_) {
  iggMarkdownLinkCallbackData dataWrapped;
  dataWrapped.link = (char *)(data_.link);
  dataWrapped.link_len = data_.linkLength;

  iggMarkdownImageData src = goMarkdownImageCallback(dataWrapped);

  ImGui::MarkdownImageData result;

  result.useLinkCallback = src.useLinkCallback;
  ImTextureID texture = reinterpret_cast<ImTextureID>(src.texture);
  if (texture == 0) {
    return result;
  }

  result.user_texture_id = texture;

  // scale image size to avoid situation, when image is larger than available region
  int availableW = ImGui::GetContentRegionAvail().x;
  if (src.shouldScale && src.size.x > availableW) {
    src.size.y = src.size.y * availableW / src.size.x;
    src.size.x = availableW;
  }

  Vec2Wrapper size(&src.size);

  Vec2Wrapper uv0(&src.uv0);
  Vec2Wrapper uv1(&src.uv1);
  Vec4Wrapper tintColor(&src.tintColor);
  Vec4Wrapper borderColor(&src.borderColor);

  result.size = *size;
  result.uv0 = *uv0;
  result.uv1 = *uv1;
  result.tint_col = *tintColor;
  result.border_col = *borderColor;

  result.isValid = true;

  return result;
}
