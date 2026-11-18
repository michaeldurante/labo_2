'''
apri webpage dove è contenuto il notebook in formato html
ctrl+shift+i per aprire ispezione pagina
cerca <iframe e tasto destro su src='link'
copy link address e incolla su nuova scheda
ctrl-s salva pagina come html
sostituisci input_file con nome del file scaricato e avvia il programma

'''
from bs4 import BeautifulSoup
import json

# HTML scaricato dall'iframe
input_file = "lecture07_differential_photometry.html"
output_file = "lectures/lecture07_differential_photometry.ipynb"

with open(input_file, "r", encoding="utf-8") as f:
    soup = BeautifulSoup(f.read(), "lxml")

notebook = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {},
    "cells": []
}

# Cerchiamo le celle del notebook
cells = soup.select("div.jp-Cell")

print(f"Celle trovate: {len(cells)}")

for d in cells:
    classes = d.get("class", [])

    # CELLE MARKDOWN
    if "jp-MarkdownCell" in classes:
        content = d.select_one(".jp-RenderedMarkdown")

        if content:
            # Rimuove i simboli ¶ dei collegamenti ai titoli
            for anchor in content.select("a.anchor-link"):
                anchor.decompose()

            notebook["cells"].append({
                "cell_type": "markdown",
                "metadata": {},
                "source": content.decode_contents()
            })

    # CELLE CODICE PYTHON
    elif "jp-CodeCell" in classes:
        editor = d.select_one(".jp-InputArea-editor")

        if editor:
            # CodeMirror: ogni riga è un elemento cm-line
            lines = editor.select(".cm-line")

            if lines:
                source = "\n".join(
                    line.get_text() for line in lines
                )
            else:
                # Strutture HTML alternative
                content = editor.select_one(".cm-content")
                if content:
                    source = content.get_text()
                else:
                    source = editor.get_text()

            notebook["cells"].append({
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": source
            })

# Salvataggio
with open(output_file, "w", encoding="utf-8") as f:
    json.dump(notebook, f, ensure_ascii=False, indent=2)

print(f"Notebook creato: {output_file}")
print(f"Totale celle convertite: {len(notebook['cells'])}")
print("Markdown:", sum(c["cell_type"] == "markdown" for c in notebook["cells"]))
print("Codice:", sum(c["cell_type"] == "code" for c in notebook["cells"]))
