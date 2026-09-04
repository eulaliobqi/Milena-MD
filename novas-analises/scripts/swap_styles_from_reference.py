#!/usr/bin/env python3
"""Pos-processamento obrigatorio depois de gerar um .docx com docx-js
(build_artigo_docx.js): substitui word/styles.xml pelo de um .docx de
referencia ja validado (ex.: docs/metodologia_resultados_DN773-GORE12T.docx).

Gotcha do docx-js: os estilos "Heading1"/"Heading2" embutidos (via
`heading: HeadingLevel.HEADING_1`) tem o `w:name` correto ("Heading 1")
mas NAO tem `<w:outlineLvl>` na definicao do estilo em styles.xml -- so
no nome. O leitor de docx do pandoc (e a navegacao/TOC do Word) usa
`w:outlineLvl` para reconhecer nivel de titulo, entao sem isso todo
paragrafo com `heading:` vira texto em negrito comum, nao um heading de
verdade (testado: nem definir `outlineLevel` direto no paragrafo, nem
sobrescrever o estilo via `Document({styles: {paragraphStyles: [...]}})`
resolve -- docx-js so duplica a definicao do estilo em vez de substituir
a default). Trocar o styles.xml inteiro por um de um .docx real
(Word/LibreOffice) e a forma confiavel de corrigir, e nao quebra nada
porque este pipeline sempre usa formatacao direta (TextRun font/size),
nunca depende do estilo para a fonte.

Uso: python scripts/swap_styles_from_reference.py <alvo.docx> <referencia.docx>
"""
import sys
import zipfile


def swap_styles(target: str, reference: str) -> None:
    with zipfile.ZipFile(reference) as z:
        ref_styles = z.read("word/styles.xml")

    with zipfile.ZipFile(target) as z:
        names = z.namelist()
        data = {n: z.read(n) for n in names}
        infos = {n: z.getinfo(n) for n in names}

    data["word/styles.xml"] = ref_styles

    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as z:
        for n in names:
            z.writestr(infos[n], data[n])

    print(f"swapped styles.xml in {target} from {reference}")


if __name__ == "__main__":
    swap_styles(sys.argv[1], sys.argv[2])
