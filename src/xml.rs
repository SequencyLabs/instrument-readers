//! Shared `roxmltree` parsing and query helpers.

use std::borrow::Cow;

use roxmltree::{Document, Node, ParsingOptions};

pub(crate) fn parse_xml(text: &str) -> Result<Document<'_>, roxmltree::Error> {
    let cleaned = text.trim_end_matches('\0').trim_start();
    Document::parse_with_options(
        cleaned,
        ParsingOptions {
            allow_dtd: true,
            ..ParsingOptions::default()
        },
    )
}

pub(crate) fn text_content<'a>(node: Node<'a, 'a>) -> Cow<'a, str> {
    let mut texts = node
        .children()
        .filter(|c| c.is_text())
        .filter_map(|c| c.text());
    let Some(first) = texts.next() else {
        return Cow::Borrowed("");
    };
    let Some(second) = texts.next() else {
        return Cow::Borrowed(first);
    };
    let mut buf = String::with_capacity(first.len() + second.len());
    buf.push_str(first);
    buf.push_str(second);
    for part in texts {
        buf.push_str(part);
    }
    Cow::Owned(buf)
}

pub(crate) fn trimmed_text(node: Node<'_, '_>) -> Option<String> {
    let text = text_content(node);
    let trimmed = text.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_owned())
}

pub(crate) fn elements<'a, 'input, 't>(
    node: Node<'a, 'input>,
    tag: &'t str,
) -> impl DoubleEndedIterator<Item = Node<'a, 'input>> + use<'a, 'input, 't> {
    node.descendants()
        .filter(move |n| n.is_element() && n.tag_name().name() == tag)
}
