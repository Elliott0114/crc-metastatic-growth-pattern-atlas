"""Manifest and workbook I/O shared by the reproduction entry point."""
from pathlib import Path
import csv
import hashlib
import json
import math
import zipfile
import xml.etree.ElementTree as ET

S = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P = "http://schemas.openxmlformats.org/package/2006/relationships"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def column(index):
    result = ""
    while index:
        index, digit = divmod(index - 1, 26)
        result = chr(65 + digit) + result
    return result


def write_workbook(path, sheets):
    """Write typed cells with the project's dependency-free OOXML approach."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        types = ET.Element("Types", xmlns="http://schemas.openxmlformats.org/package/2006/content-types")
        ET.SubElement(types, "Default", Extension="rels", ContentType="application/vnd.openxmlformats-package.relationships+xml")
        ET.SubElement(types, "Default", Extension="xml", ContentType="application/xml")
        ET.SubElement(types, "Override", PartName="/xl/workbook.xml", ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml")
        book = ET.Element("workbook", xmlns=S, attrib={"xmlns:r": R})
        entries = ET.SubElement(book, "sheets")
        links = ET.Element("Relationships", xmlns=P)
        root = ET.Element("Relationships", xmlns=P)
        ET.SubElement(root, "Relationship", Id="rId1", Type=R + "/officeDocument", Target="xl/workbook.xml")
        archive.writestr("_rels/.rels", ET.tostring(root))
        for index, (name, rows) in enumerate(sheets, 1):
            assert len(name) <= 31
            ET.SubElement(entries, "sheet", name=name, sheetId=str(index), attrib={"r:id": f"rId{index}"})
            target = f"worksheets/sheet{index}.xml"
            ET.SubElement(links, "Relationship", Id=f"rId{index}", Type=R + "/worksheet", Target=target)
            ET.SubElement(types, "Override", PartName="/xl/" + target, ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml")
            sheet = ET.Element("worksheet", xmlns=S)
            view = ET.SubElement(ET.SubElement(sheet, "sheetViews"), "sheetView", workbookViewId="0")
            ET.SubElement(view, "pane", ySplit="1", topLeftCell="A2", activePane="bottomLeft", state="frozen")
            cols = ET.SubElement(sheet, "cols")
            ET.SubElement(cols, "col", min="1", max=str(max(map(len, rows), default=1)), width="24", customWidth="1")
            data = ET.SubElement(sheet, "sheetData")
            for ri, values in enumerate(rows, 1):
                row = ET.SubElement(data, "row", r=str(ri))
                for ci, value in enumerate(values, 1):
                    cell = ET.SubElement(row, "c", r=f"{column(ci)}{ri}")
                    if value is None:
                        continue
                    if isinstance(value, bool):
                        cell.set("t", "b")
                        ET.SubElement(cell, "v").text = "1" if value else "0"
                    elif isinstance(value, (int, float)) and math.isfinite(value):
                        ET.SubElement(cell, "v").text = format(value, ".17g")
                    else:
                        cell.set("t", "inlineStr")
                        text = ET.SubElement(ET.SubElement(cell, "is"), "t")
                        text.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
                        text.text = str(value)
            archive.writestr("xl/" + target, ET.tostring(sheet, encoding="utf-8"))
        archive.writestr("[Content_Types].xml", ET.tostring(types))
        archive.writestr("xl/workbook.xml", ET.tostring(book, encoding="utf-8"))
        archive.writestr("xl/_rels/workbook.xml.rels", ET.tostring(links))


def read_workbook(path):
    result = []
    ns = {"m": S}
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = ["".join(x.itertext()) for x in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        links = {x.attrib["Id"]: x.attrib["Target"] for x in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))}
        for sheet in ET.fromstring(archive.read("xl/workbook.xml")).findall("m:sheets/m:sheet", ns):
            target = links[sheet.attrib["{" + R + "}id"]]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            rows = []
            for row in ET.fromstring(archive.read(target)).findall("m:sheetData/m:row", ns):
                values = []
                for cell in row.findall("m:c", ns):
                    col = 0
                    for letter in (x for x in cell.attrib["r"] if x.isalpha()):
                        col = col * 26 + ord(letter) - 64
                    values.extend([None] * (col - len(values)))
                    typ = cell.attrib.get("t")
                    value = cell.find("m:v", ns)
                    if typ == "inlineStr":
                        value = "".join(x.text or "" for x in cell.findall("m:is//m:t", ns))
                    elif value is None:
                        value = None
                    elif typ == "s":
                        value = strings[int(value.text)]
                    elif typ == "b":
                        value = value.text == "1"
                    elif typ in ("str", "e"):
                        value = value.text
                    else:
                        value = float(value.text)
                    values[col - 1] = value
                rows.append(values)
            result.append((sheet.attrib["name"], rows))
    return result
