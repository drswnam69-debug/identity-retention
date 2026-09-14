"""Assemble the bioRxiv preprint PDF and a DOCX from the v2 master markdown."""
import re
import subprocess
import pathlib

SRC = pathlib.Path("manuscript_v4.md")
BUILD = pathlib.Path("build")
BUILD.mkdir(exist_ok=True)

s = SRC.read_text()

# ---------------------------------------------------------------- front matter
# Pandoc gets title/author from the markdown body itself; strip the manual
# heading block and rebuild it as plain paragraphs so it typesets as a title page.


# ---------------------------------------------------------------- figure images
# Place each figure immediately after its legend paragraph.
for n in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12):
    pat = re.compile(r"(\*\*Figure " + str(n) + r"\..*?)(?=\n\n)", re.S)
    m = pat.search(s)
    assert m, f"legend for Figure {n} not found"
    img = f"\n\n![](build/Figure{n}.png){{width=98%}}\n"
    s = s[:m.end()] + img + s[m.end():]

# ---------------------------------------------------------------- horizontal rules
# The markdown uses --- as section separators; in LaTeX these become page-wide
# rules that add nothing. Drop them.
s = re.sub(r"\n---\n", "\n", s)

PRE = BUILD / "preprint_body.md"
PRE.write_text(s)

HEADER = BUILD / "header.tex"
HEADER.write_text(r"""
\usepackage{longtable,booktabs}
\usepackage{float}
\makeatletter
\def\fps@figure{H}
\makeatother
\usepackage{ragged2e}
\AtBeginEnvironment{longtable}{\footnotesize}
\setlength{\emergencystretch}{3em}
\usepackage{etoolbox}
\preto{\longtable}{\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt}}
""")

common = [
    "pandoc", str(PRE),
    "--from", "markdown+pipe_tables+tex_math_dollars",
    "--resource-path", ".",
    "-V", "geometry:margin=2.4cm",
    "-V", "fontsize=11pt",
    "-V", "linestretch=1.25",
    "-V", "colorlinks=true",
    "-V", "linkcolor=black",
    "-V", "urlcolor=black",
]

pdf = common + [
    "--pdf-engine", "xelatex",
    "-V", "mainfont=DejaVu Serif",
    "-V", "sansfont=DejaVu Sans",
    "-V", "monofont=DejaVu Sans Mono",
    "-H", str(HEADER),
    "-o", "build/Nam_2026_preprint_v4.pdf",
]
docx = common[:4] + ["--reference-doc", "build/reference.docx",
        "-o", "build/Nam_2026_manuscript_v4.docx"]

for cmd in (pdf, docx):
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(cmd[-1], "->", r.returncode)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])

# ---------------------------------------------------------------- docx polish
# GigaScience asks for double line spacing, page numbers and continuous line
# numbers in the submitted manuscript file. Pandoc's writer emits none of the
# three, so patch them into the package directly.
import zipfile, shutil, os

FOOTER_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
    '<w:p><w:pPr><w:jc w:val="center"/><w:spacing w:after="0" w:line="240" '
    'w:lineRule="auto"/></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
    '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate"/></w:r>'
    '<w:r><w:t>1</w:t></w:r>'
    '<w:r><w:fldChar w:fldCharType="end"/></w:r>'
    '</w:p></w:ftr>'
)

FOOTER_RID = "rId9001"
SECTPR = (
    '<w:sectPr>'
    '<w:footerReference w:type="default" r:id="%s"/>'
    '<w:pgSz w:w="11906" w:h="16838"/>'
    '<w:pgMar w:top="1418" w:right="1418" w:bottom="1418" w:left="1418" '
    'w:header="709" w:footer="709" w:gutter="0"/>'
    '<w:lnNumType w:countBy="1" w:restart="continuous" w:distance="360"/>'
    '</w:sectPr>' % FOOTER_RID
)


def polish_docx(path):
    tmp = path + "._x"
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)
    with zipfile.ZipFile(path) as z:
        z.extractall(tmp)

    # -- double spacing -------------------------------------------------
    sp = os.path.join(tmp, "word", "styles.xml")
    t = open(sp, encoding="utf-8").read()
    t = t.replace('<w:spacing w:after="200" />',
                  '<w:spacing w:after="200" w:line="480" w:lineRule="auto" />')
    open(sp, "w", encoding="utf-8").write(t)

    # -- footer part ----------------------------------------------------
    open(os.path.join(tmp, "word", "footer1.xml"), "w", encoding="utf-8").write(FOOTER_XML)

    ctp = os.path.join(tmp, "[Content_Types].xml")
    ct = open(ctp, encoding="utf-8").read()
    assert "footer1.xml" not in ct
    ct = ct.replace(
        "</Types>",
        '<Override PartName="/word/footer1.xml" ContentType="application/vnd.'
        'openxmlformats-officedocument.wordprocessingml.footer+xml" /></Types>')
    open(ctp, "w", encoding="utf-8").write(ct)

    rp = os.path.join(tmp, "word", "_rels", "document.xml.rels")
    rels = open(rp, encoding="utf-8").read()
    assert FOOTER_RID not in rels
    rels = rels.replace(
        "</Relationships>",
        '<Relationship Type="http://schemas.openxmlformats.org/officeDocument/'
        '2006/relationships/footer" Id="%s" Target="footer1.xml" />'
        "</Relationships>" % FOOTER_RID)
    open(rp, "w", encoding="utf-8").write(rels)

    # -- section properties: page size, margins, line numbers, footer ---
    dp = os.path.join(tmp, "word", "document.xml")
    doc = open(dp, encoding="utf-8").read()
    n = doc.count("<w:sectPr />")
    assert n == 1, "expected exactly one empty sectPr, found %d" % n
    doc = doc.replace("<w:sectPr />", SECTPR)
    open(dp, "w", encoding="utf-8").write(doc)

    # -- repack ---------------------------------------------------------
    os.remove(path)
    zf = zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(tmp):
        for f in files:
            full = os.path.join(root, f)
            zf.write(full, os.path.relpath(full, tmp))
    zf.close()
    shutil.rmtree(tmp)
    print("polished (double spacing, page numbers, line numbers):", path)


polish_docx("build/Nam_2026_manuscript_v4.docx")
