#!/usr/bin/env python3
"""
OpenWebUI plugin code audit.
Verifies that tool/function Python code stored in webui.db compiles and imports cleanly.
"""

import argparse
import ast
import importlib.util
import sqlite3
import sys


def find_missing_imports(source: str, import_check: bool = True) -> list:
    """Find modules that are imported but not available."""
    if not import_check:
        return []

    missing = []
    seen = set()
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return missing

    for node in ast.walk(tree):
        module_name = ""
        line_no = getattr(node, "lineno", 0) or 0

        if isinstance(node, ast.Import):
            for imported in node.names:
                module_name = (imported.name or "").split(".")[0]
                if not module_name or module_name in seen:
                    continue
                if importlib.util.find_spec(module_name) is None:
                    missing.append((module_name, line_no))
                    seen.add(module_name)

        if isinstance(node, ast.ImportFrom):
            if node.level and node.level > 0:
                continue
            module_name = (node.module or "").split(".")[0]
            if not module_name or module_name in seen:
                continue
            if importlib.util.find_spec(module_name) is None:
                missing.append((module_name, line_no))
                seen.add(module_name)

    return missing


def audit_plugins(focus: str, import_check: bool, db_path: str = "/app/backend/data/webui.db") -> tuple:
    """
    Audit plugin code from webui.db.

    Returns:
        tuple: (checked_count, issues_list)
    """
    tables = []
    if focus in ("all", "tool"):
        tables.append("tool")
    if focus in ("all", "function"):
        tables.append("function")

    con = sqlite3.connect(db_path)
    cur = con.cursor()

    checked = 0
    issues = []

    for table in tables:
        try:
            cur.execute(f"SELECT id, name, content FROM {table}")
        except sqlite3.OperationalError:
            continue

        for row_id, name, content in cur.fetchall():
            checked += 1
            source = content or ""
            try:
                compile(source, f"{table}:{row_id}", "exec")
            except SyntaxError as exc:
                lines = source.splitlines()
                line_txt = ""
                if exc.lineno and 1 <= exc.lineno <= len(lines):
                    line_txt = lines[exc.lineno - 1]
                issues.append(
                    {
                        "table": table,
                        "id": row_id,
                        "name": name,
                        "line": exc.lineno or 0,
                        "offset": exc.offset or 0,
                        "msg": str(exc.msg),
                        "code": line_txt,
                    }
                )
                continue

            for module_name, line_no in find_missing_imports(source, import_check):
                issues.append(
                    {
                        "table": table,
                        "id": row_id,
                        "name": name,
                        "line": line_no,
                        "offset": 0,
                        "msg": f"missing import dependency: {module_name}",
                        "code": "",
                    }
                )

    con.close()
    return checked, issues


def main():
    parser = argparse.ArgumentParser(description="Audit OpenWebUI plugin code")
    parser.add_argument("focus", choices=["all", "tool", "function"], help="Scope of audit")
    parser.add_argument("import_check", nargs="?", default="true", help="Check import dependencies (true/false)")
    parser.add_argument("--db-path", default="/app/backend/data/webui.db", help="Path to webui.db")
    args = parser.parse_args()

    import_check = args.import_check.lower() in ("1", "true", "yes")
    checked, issues = audit_plugins(args.focus, import_check, args.db_path)

    print(f"AUDIT_CHECKED={checked}")
    print(f"AUDIT_ISSUES={len(issues)}")
    for issue in issues:
        print(
            "ISSUE|{table}|{id}|{name}|{line}|{offset}|{msg}|{code}".format(
                table=issue["table"],
                id=issue["id"],
                name=(issue["name"] or "").replace("\n", " "),
                line=issue["line"],
                offset=issue["offset"],
                msg=issue["msg"].replace("\n", " "),
                code=issue["code"].replace("\n", " "),
            )
        )

    sys.exit(1 if issues else 0)


if __name__ == "__main__":
    main()
