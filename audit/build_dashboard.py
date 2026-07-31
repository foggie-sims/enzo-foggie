#!/usr/bin/env python3
"""Render audit/dashboard.html from audit/status.yml.

Usage:
    python3 audit/build_dashboard.py            # local rebuild
    CI_RUN_URL=... CI_RUN_STATUS=... CI_RUN_SHA=... CI_RUN_TIME=... \
        python3 audit/build_dashboard.py        # CI rebuild with run info

Requires PyYAML. Output is pure-ASCII HTML (non-ASCII escaped to entities) so
it renders correctly under any charset, viewed locally or hosted.
"""

import datetime
import html
import os
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")

HERE = os.path.dirname(os.path.abspath(__file__))
STATUS_FILE = os.path.join(HERE, "status.yml")
OUT_FILE = os.path.join(HERE, "dashboard.html")

REPO_URL = "https://github.com/foggie-sims/enzo-foggie"
AUDIT_FILE = "enzoaudit.html"  # sits next to the dashboard in audit/

STATUS_ORDER = ["done", "in-progress", "todo", "wontfix"]
STATUS_LABEL = {"done": "Done", "in-progress": "In progress",
                "todo": "To do", "wontfix": "Won't fix"}


def esc(s):
    """HTML-escape and force ASCII with numeric entities."""
    return html.escape(str(s)).encode("ascii", "xmlcharrefreplace").decode("ascii")


def ref_anchor(ref):
    """Map an audit section ref like '3.2' to the anchor id in enzoaudit.html."""
    ref = str(ref)
    if ref.upper().startswith("A"):
        return "appendix"
    return "s" + ref.replace(".", "")


def render_item(item):
    status = item.get("status", "todo")
    chip = f'<span class="chip st-{esc(status)}">{esc(STATUS_LABEL.get(status, status))}</span>'
    klass = item.get("class", "")
    klass_chip = f'<span class="chip cls-{esc(klass).lower()}">{esc(klass)}</span>' if klass else ""
    ref = item.get("ref", "")
    ref_link = (f'<a href="{AUDIT_FILE}#{ref_anchor(ref)}">&#167;{esc(ref)}</a>'
                if ref else "&#8212;")
    prs = item.get("prs") or []
    pr_links = ", ".join(
        f'<a href="{REPO_URL}/pull/{int(n)}">#{int(n)}</a>' for n in prs) or "&#8212;"
    tests = item.get("tests") or []
    test_bits = []
    for t in tests:
        t = str(t)
        if t.startswith("http://") or t.startswith("https://"):
            test_bits.append(f'<a href="{esc(t)}">{esc(t.split("/")[-1] or "link")}</a>')
        else:
            test_bits.append(esc(t))
    tests_html = "; ".join(test_bits) or "&#8212;"
    note = item.get("note", "")
    note_html = f'<div class="note-line">{esc(note)}</div>' if note else ""
    return (f"<tr class=\"row-{esc(status)}\">"
            f"<td>{esc(item.get('id', ''))}</td>"
            f"<td>{esc(item.get('title', ''))}{note_html}</td>"
            f"<td>{ref_link}</td>"
            f"<td>{klass_chip}</td>"
            f"<td>{chip}</td>"
            f"<td>{pr_links}</td>"
            f"<td>{tests_html}</td>"
            f"</tr>")


def tier_counts(items):
    c = {k: 0 for k in STATUS_ORDER}
    for it in items:
        c[it.get("status", "todo")] = c.get(it.get("status", "todo"), 0) + 1
    return c


def progress_bar(counts, total):
    if total == 0:
        return ""
    done = counts.get("done", 0)
    prog = counts.get("in-progress", 0)
    dw = 100.0 * done / total
    pw = 100.0 * prog / total
    return (f'<div class="bar" role="img" aria-label="{done} of {total} done">'
            f'<div class="bar-done" style="width:{dw:.1f}%"></div>'
            f'<div class="bar-prog" style="width:{pw:.1f}%"></div></div>')


def main():
    with open(STATUS_FILE, encoding="utf-8") as f:
        data = yaml.safe_load(f)

    tiers = data.get("tiers", [])
    all_items = [it for t in tiers for it in t.get("items", [])]
    totals = tier_counts(all_items)
    n_all = len(all_items)

    gen_time = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    ci_url = os.environ.get("CI_RUN_URL", "")
    ci_status = os.environ.get("CI_RUN_STATUS", "")
    ci_sha = os.environ.get("CI_RUN_SHA", "")[:9]
    ci_time = os.environ.get("CI_RUN_TIME", "")
    if ci_url:
        badge_class = {"success": "ok", "failure": "bad"}.get(ci_status, "meh")
        ci_block = (f'<p class="ci-line">Latest branch CI: '
                    f'<a class="ci-badge {badge_class}" href="{esc(ci_url)}">{esc(ci_status or "unknown")}</a>'
                    f' at <code>{esc(ci_sha)}</code> ({esc(ci_time)}) &#183; '
                    f'<a href="{REPO_URL}/actions?query=branch%3Aenzo-performance">all runs</a></p>')
    else:
        ci_block = (f'<p class="ci-line">Generated locally (no CI context) &#183; '
                    f'<a href="{REPO_URL}/actions?query=branch%3Aenzo-performance">branch CI runs</a></p>')

    # Headline measured results (status.yml: findings). These are the
    # conclusions the ledger rows cannot convey on their own - several
    # overturn the original audit's ranking, so they lead the page.
    findings = data.get("findings") or {}
    fitems = findings.get("items") or []
    if fitems:
        cards = "\n".join(
            f'<div class="finding"><div class="finding-title">{esc(fi.get("title",""))}</div>'
            f'<div class="finding-detail">{esc(fi.get("detail",""))}</div></div>'
            for fi in fitems)
        anchor = findings.get("anchor", "")
        findings_block = (
            '<section class="findings">'
            '<h2 id="findings">Measured findings</h2>'
            + (f'<p class="tier-desc">Anchor: {esc(anchor)}</p>' if anchor else "")
            + f'<div class="finding-grid">{cards}</div>'
            '<p class="tier-desc">Full prioritised plan: '
            '<a href="ROADMAP.md">ROADMAP.md</a> &#183; '
            'deferred barrier work: <a href="SPFinalize_Edits.md">SPFinalize_Edits.md</a></p>'
            '</section>')
    else:
        findings_block = ""

    tier_sections = []
    summary_cards = []
    for t in tiers:
        items = t.get("items", [])
        counts = tier_counts(items)
        total = len(items)
        done = counts.get("done", 0)
        summary_cards.append(
            f'<a class="card" href="#{esc(t["id"])}">'
            f'<div class="card-name">{esc(t["name"])}</div>'
            f'<div class="card-num"><span class="big">{done}</span>/{total} done'
            + (f' &#183; {counts["in-progress"]} in progress' if counts.get("in-progress") else "")
            + f'</div>{progress_bar(counts, total)}</a>')
        rows = "\n".join(render_item(it) for it in items)
        tier_sections.append(f"""
  <h2 id="{esc(t['id'])}">{esc(t['name'])}</h2>
  <p class="tier-desc">{esc(t.get('description', ''))} &#8212; {done}/{total} done</p>
  <div class="tablewrap">
    <table>
      <thead><tr><th>ID</th><th>Item</th><th>Audit</th><th>Class</th><th>Status</th><th>PRs</th><th>Test evidence</th></tr></thead>
      <tbody>
{rows}
      </tbody>
    </table>
  </div>""")

    page = f"""<title>Enzo/FOGGIE Audit Remediation Dashboard</title>
<style>
  :root {{
    --paper:#faf9f7; --ink:#1c1b22; --mist:#5f5c6b; --violet:#54468a;
    --violet-soft:#eceafa; --amber:#a34d09; --amber-soft:#f7ede2;
    --line:#e3e1da; --card:#f3f2ee; --code-bg:#edecf2;
    --good:#2c6e49; --good-soft:#e4f0e9; --bad:#a3341e;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --paper:#16151c; --ink:#e8e6f0; --mist:#a09db0; --violet:#a99ee2;
      --violet-soft:#262336; --amber:#e0a35c; --amber-soft:#2e2418;
      --line:#2d2b38; --card:#1e1d26; --code-bg:#24222e;
      --good:#7fc8a0; --good-soft:#1c2b22; --bad:#e08a76; }}
  }}
  :root[data-theme="dark"] {{ --paper:#16151c; --ink:#e8e6f0; --mist:#a09db0; --violet:#a99ee2;
    --violet-soft:#262336; --amber:#e0a35c; --amber-soft:#2e2418;
    --line:#2d2b38; --card:#1e1d26; --code-bg:#24222e;
    --good:#7fc8a0; --good-soft:#1c2b22; --bad:#e08a76; }}
  :root[data-theme="light"] {{ --paper:#faf9f7; --ink:#1c1b22; --mist:#5f5c6b; --violet:#54468a;
    --violet-soft:#eceafa; --amber:#a34d09; --amber-soft:#f7ede2;
    --line:#e3e1da; --card:#f3f2ee; --code-bg:#edecf2;
    --good:#2c6e49; --good-soft:#e4f0e9; --bad:#a3341e; }}
  html {{ background: var(--paper); }}
  body {{ margin:0; background:var(--paper); color:var(--ink);
    font-family: Seravek, Avenir, "Avenir Next", "Segoe UI", system-ui, sans-serif;
    font-size:0.95rem; line-height:1.55; }}
  .wrap {{ max-width: 62rem; margin:0 auto; padding: 2.5rem 1.25rem 5rem; }}
  code {{ font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
    font-size:0.85em; background:var(--code-bg); border-radius:3px; padding:0.08em 0.35em; }}
  .eyebrow {{ text-transform:uppercase; letter-spacing:0.14em; font-size:0.7rem;
    font-weight:600; color:var(--amber); margin:0 0 0.6rem; }}
  h1 {{ font-size: clamp(1.5rem, 4vw, 2.1rem); margin:0 0 0.4rem; font-weight:700; }}
  .gen-line, .ci-line {{ color:var(--mist); font-size:0.8rem; margin:0.15rem 0; }}
  .ci-badge {{ font-weight:700; text-transform:uppercase; font-size:0.7rem;
    letter-spacing:0.06em; border-radius:999px; padding:0.12rem 0.6rem; text-decoration:none; }}
  .ci-badge.ok  {{ background:var(--good); color:var(--paper); }}
  .ci-badge.bad {{ background:var(--bad); color:var(--paper); }}
  .ci-badge.meh {{ background:var(--mist); color:var(--paper); }}
  .cards {{ display:grid; grid-template-columns:repeat(auto-fit, minmax(13rem,1fr));
    gap:0.8rem; margin:1.6rem 0 1rem; }}
  .card {{ background:var(--card); border-radius:8px; padding:0.9rem 1rem;
    text-decoration:none; color:var(--ink); }}
  .card:hover {{ outline:2px solid var(--violet); }}
  .card-name {{ font-weight:650; font-size:0.82rem; margin-bottom:0.35rem; }}
  .card-num {{ color:var(--mist); font-size:0.8rem; margin-bottom:0.5rem; }}
  .card-num .big {{ color:var(--ink); font-weight:700; font-size:1.25rem; }}
  .bar {{ display:flex; height:0.45rem; border-radius:999px; overflow:hidden;
    background:var(--line); }}
  .bar-done {{ background:var(--good); }}
  .bar-prog {{ background:var(--amber); }}
  h2 {{ font-size:1.2rem; margin:2.6rem 0 0.3rem; padding-top:1.1rem;
    border-top:2px solid var(--violet); }}
  .tier-desc {{ color:var(--mist); font-size:0.85rem; margin:0 0 0.8rem; }}
  .tablewrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:8px; }}
  table {{ border-collapse:collapse; width:100%; font-size:0.83rem; }}
  th {{ text-align:left; font-size:0.68rem; text-transform:uppercase; letter-spacing:0.07em;
    color:var(--mist); background:var(--card); padding:0.5rem 0.7rem;
    border-bottom:1px solid var(--line); white-space:nowrap; }}
  td {{ padding:0.5rem 0.7rem; border-bottom:1px solid var(--line); vertical-align:top; }}
  tr:last-child td {{ border-bottom:none; }}
  td:first-child {{ font-weight:650; white-space:nowrap; }}
  tr.row-done td {{ opacity:0.75; }}
  .note-line {{ color:var(--mist); font-size:0.78rem; margin-top:0.2rem; }}
  .findings {{ margin:1.6rem 0; }}
  .finding-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(280px,1fr));
                   gap:0.75rem; margin-top:0.6rem; }}
  .finding {{ border:1px solid var(--line); border-radius:8px; padding:0.7rem 0.85rem;
              background:var(--card); }}
  .finding-title {{ font-weight:600; font-size:0.9rem; margin-bottom:0.3rem; }}
  .finding-detail {{ color:var(--mist); font-size:0.8rem; line-height:1.45; }}
  .chip {{ display:inline-block; font-size:0.64rem; font-weight:700; text-transform:uppercase;
    letter-spacing:0.07em; border-radius:999px; padding:0.1rem 0.55rem; white-space:nowrap; }}
  .chip.st-done {{ background:var(--good); color:var(--paper); }}
  .chip.st-in-progress {{ background:var(--amber); color:var(--paper); }}
  .chip.st-todo {{ background:transparent; color:var(--mist); border:1px solid var(--mist); }}
  .chip.st-wontfix {{ background:var(--mist); color:var(--paper); }}
  .chip.cls-c0 {{ background:var(--good-soft); color:var(--good); }}
  .chip.cls-c1 {{ background:var(--violet-soft); color:var(--violet); }}
  .chip.cls-c2 {{ background:var(--amber-soft); color:var(--amber); }}
  a {{ color:var(--violet); }}
  a:hover {{ color:var(--amber); }}
  footer {{ margin-top:3rem; border-top:1px solid var(--line); padding-top:1rem;
    font-size:0.78rem; color:var(--mist); }}
</style>
<div class="wrap">
  <p class="eyebrow">enzo-performance &#183; audit remediation</p>
  <h1>Audit Remediation Dashboard</h1>
  <p class="gen-line">Generated {esc(gen_time)} from <code>audit/status.yml</code>
    &#183; overall {totals.get('done', 0)}/{n_all} done, {totals.get('in-progress', 0)} in progress</p>
  {ci_block}
  <div class="cards">
{''.join(summary_cards)}
  </div>
{findings_block}
{''.join(tier_sections)}
  <footer>
    <p>Rebuilt automatically on every change to <code>audit/status.yml</code> on the
    <code>enzo-performance</code> branch. To claim an item, edit the ledger in the PR
    that fixes it: set its status, add the PR number, and link test evidence.
    Item text is abbreviated &#8212; the full recommendation is in
    <a href="{AUDIT_FILE}">the audit</a> (&#167; column) and the testing standard per
    class (C0/C1/C2) is in <a href="enzotestplan.html">the test plan</a>.</p>
  </footer>
</div>
"""
    with open(OUT_FILE, "w", encoding="ascii") as f:
        f.write(page)
    print(f"Wrote {OUT_FILE}: {n_all} items, {totals.get('done', 0)} done, "
          f"{totals.get('in-progress', 0)} in progress")


if __name__ == "__main__":
    main()
