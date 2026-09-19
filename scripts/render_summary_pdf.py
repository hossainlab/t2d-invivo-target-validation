import os
import subprocess
import markdown

def main():
    md_path = os.path.abspath("figures/Executive_Summary.md")
    with open(md_path, "r", encoding="utf-8") as f:
        md_content = f.read()

    html_body = markdown.markdown(md_content, extensions=["tables", "fenced_code", "nl2br"])

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Executive Summary & Central Conclusion</title>
<style>
  @page {
    size: A4 portrait;
    margin: 12mm 14mm 12mm 14mm;
  }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    color: #0f172a;
    line-height: 1.4;
    font-size: 9.5pt;
    background-color: #ffffff;
    margin: 0;
    padding: 0;
  }
  h1 {
    font-size: 16pt;
    font-weight: 800;
    color: #0f172a;
    border-bottom: 2px solid #2563eb;
    padding-bottom: 4px;
    margin-top: 0;
    margin-bottom: 4px;
    letter-spacing: -0.02em;
  }
  h2 {
    font-size: 11.5pt;
    font-weight: 750;
    color: #1e293b;
    border-bottom: 1.2px solid #cbd5e1;
    padding-bottom: 2px;
    margin-top: 12px;
    margin-bottom: 6px;
    page-break-after: avoid;
  }
  h3 {
    font-size: 10pt;
    font-weight: 650;
    color: #334155;
    margin-top: 8px;
    margin-bottom: 3px;
    page-break-after: avoid;
  }
  p {
    margin-top: 2px;
    margin-bottom: 5px;
  }
  strong {
    color: #0f172a;
  }
  hr {
    border: 0;
    height: 1px;
    background: #e2e8f0;
    margin: 8px 0;
  }
  pre {
    background-color: #f8fafc;
    border: 1px solid #cbd5e1;
    border-radius: 4px;
    padding: 6px 8px;
    font-family: "Cascadia Code", "Consolas", "Courier New", monospace;
    font-size: 7.2pt;
    line-height: 1.2;
    overflow-x: auto;
    color: #0f172a;
    margin: 6px 0;
    page-break-inside: avoid;
  }
  code {
    font-family: "Cascadia Code", "Consolas", monospace;
    background-color: #f1f5f9;
    padding: 1px 3px;
    border-radius: 3px;
    font-size: 8.5pt;
    color: #0369a1;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    margin: 6px 0;
    font-size: 8pt;
    page-break-inside: avoid;
  }
  th, td {
    border: 1px solid #cbd5e1;
    padding: 4.5px 6px;
    text-align: left;
    vertical-align: top;
  }
  th {
    background-color: #f1f5f9;
    font-weight: 650;
    color: #1e293b;
  }
  tr:nth-child(even) td {
    background-color: #f8fafc;
  }
  blockquote {
    border-left: 3.5px solid #2563eb;
    background-color: #eff6ff;
    padding: 6px 10px;
    margin: 6px 0;
    border-radius: 0 4px 4px 0;
    font-size: 9pt;
    page-break-inside: avoid;
  }
  ul, ol {
    margin-top: 2px;
    margin-bottom: 4px;
    padding-left: 18px;
  }
  li {
    margin-bottom: 2px;
  }
</style>
</head>
<body>
""" + html_body + """
</body>
</html>"""

    html_file = os.path.abspath("figures/Executive_Summary.html")
    with open(html_file, "w", encoding="utf-8") as f:
        f.write(html_template)
    print(f"Wrote {html_file}")

    chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
    edge_path = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
    browser = chrome_path if os.path.exists(chrome_path) else edge_path

    out_pdf = os.path.abspath("figures/Executive_Summary.pdf")
    cmd = [
        browser,
        "--headless",
        "--disable-gpu",
        "--allow-file-access-from-files",
        "--no-pdf-header-footer",
        f"--print-to-pdf={out_pdf}",
        html_file
    ]

    print("Rendering PDF...")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0 and os.path.exists(out_pdf):
        print(f"Successfully generated {out_pdf} ({os.path.getsize(out_pdf)} bytes)")
    else:
        print(f"Error rendering PDF: {res.stderr}")

if __name__ == "__main__":
    main()
