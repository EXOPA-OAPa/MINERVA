"""
Check that two versions of a Python file contain the same code, ignoring
comments, docstrings, `pass` statements next to other
statements and human-readable messages.

Strings are ignored only where they are messages: inside print(...), inside
exception constructors in `raise`, in assignments to variables named
`*_msg` (and harm_str, status_str), in list.append on lists named `*_reasons`, and text arguments of matplotlib calls
(set_title, set_xlabel, set_ylabel, suptitle, text, annotate, legend, and the
`label=`/`title=` keywords). Every other string (dictionary keys, tags such as
'slow', fit descriptions compared in the code, file names...) must be identical.

Usage:  python scripts/check_same_code.py old.py new.py
"""
import ast
import sys

TEXT_METHODS = {'set_title', 'set_xlabel', 'set_ylabel', 'suptitle', 'text', 'annotate',
                'legend', 'title', 'xlabel', 'ylabel', 'set_label'}
TEXT_KEYWORDS = {'label', 'title', 'desc'}   # desc: progress-bar text (tqdm)
# variables that only hold log messages (besides the *_msg ones)
MESSAGE_VARS = {'harm_str', 'status_str'}


class Normalise(ast.NodeTransformer):
    def _blank_strings(self, node):
        for sub in ast.walk(node):
            if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                sub.value = '<msg>'
            if isinstance(sub, ast.JoinedStr):
                sub.values = [ast.Constant(value='<msg>')]
        return node

    def visit_Call(self, node):
        self.generic_visit(node)
        f = node.func
        name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else None
        # messages collected in lists such as rejection_reasons.append(...)
        if name == 'append' and isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) \
                and f.value.id.endswith('_reasons'):
            for a in node.args:
                self._blank_strings(a)
        # screen messages: sys.stdout.write(...) / sys.__stdout__.write(...) (never file writes)
        if name == 'write' and isinstance(f, ast.Attribute) and isinstance(f.value, ast.Attribute) \
                and isinstance(f.value.value, ast.Name) and f.value.value.id == 'sys' \
                and f.value.attr in ('stdout', '__stdout__'):
            for a in node.args:
                self._blank_strings(a)
        if name == 'print' or name in TEXT_METHODS:
            for a in node.args:
                self._blank_strings(a)
        for kw in node.keywords:
            if kw.arg in TEXT_KEYWORDS:
                self._blank_strings(kw.value)
        return node

    def visit_Assign(self, node):
        self.generic_visit(node)
        # variables named *_msg hold messages (printed and/or written to the log)
        if all(isinstance(t, ast.Name) and (t.id.endswith('_msg') or t.id in MESSAGE_VARS) for t in node.targets):
            self._blank_strings(node.value)
        return node

    def visit_Raise(self, node):
        self.generic_visit(node)
        if node.exc is not None:
            self._blank_strings(node.exc)
        return node

    def _drop_docstring(self, node):
        self.generic_visit(node)
        if node.body and isinstance(node.body[0], ast.Expr) and isinstance(getattr(node.body[0], 'value', None), ast.Constant) \
                and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:] or [ast.Pass()]
        return node

    def generic_visit(self, node):
        super().generic_visit(node)
        # a `pass` next to other statements, or a bare string (docstring or disabled code
        # in triple quotes), does nothing: ignore it
        for field in ('body', 'orelse', 'finalbody'):
            stmts = getattr(node, field, None)
            if isinstance(stmts, list) and stmts and isinstance(stmts[0], ast.stmt):
                kept = [s for s in stmts if not isinstance(s, ast.Pass)
                        and not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant) and isinstance(s.value.value, str))]
                setattr(node, field, kept or [ast.Pass()])
        return node

    visit_FunctionDef = _drop_docstring
    visit_Module = _drop_docstring
    visit_ClassDef = _drop_docstring


def normalised_dump(path):
    tree = ast.parse(open(path, encoding='utf-8').read())
    return ast.dump(Normalise().visit(tree), include_attributes=False)


def main(old, new):
    a, b = normalised_dump(old), normalised_dump(new)
    if a == b:
        print(f'SAME CODE: {new}')
        return 0
    # report the first difference
    i = next(k for k in range(min(len(a), len(b))) if a[k] != b[k])
    print(f'CODE DIFFERS: {new}\n  old: ...{a[max(0, i - 120):i + 120]}...\n  new: ...{b[max(0, i - 120):i + 120]}...')
    return 1


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
