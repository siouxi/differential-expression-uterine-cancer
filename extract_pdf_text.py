from pathlib import Path

from pypdf import PdfReader


PDF_PATHS = [
    Path("ESTUDIOS/PAPERS/Lectura1-TMM.pdf"),
    Path("ESTUDIOS/PAPERS/Lectura2-DEGvisualization-RNAseq.pdf"),
    Path("ESTUDIOS/PAPERS/Lectura5-Enrichment.pdf"),
    Path("ESTUDIOS/PRESENTACIONES/2aIntroBioconductor-nuevo.pdf"),
    Path("ESTUDIOS/PRESENTACIONES/2bRNAseq-counttable.pdf"),
    Path("ESTUDIOS/PRESENTACIONES/3ExpresionDiferencial-nuevo.pdf"),
    Path("ESTUDIOS/PRESENTACIONES/4GenomicCategorical-nuevo.pdf"),
]


def extract_pdf_text(pdf_path: Path) -> str:
    reader = PdfReader(str(pdf_path))
    return "\n".join((page.extract_text() or "").strip() for page in reader.pages)


def main() -> None:
    for pdf_path in PDF_PATHS:
        if not pdf_path.exists():
            print(f"Skipping missing PDF: {pdf_path}")
            continue

        text = extract_pdf_text(pdf_path)
        out_path = pdf_path.with_suffix(pdf_path.suffix + ".txt")
        out_path.write_text(text, encoding="utf-8")
        print(f"Wrote text to {out_path}")


if __name__ == "__main__":
    main()

