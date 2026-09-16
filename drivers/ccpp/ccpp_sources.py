from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def sources():
    candidates = sorted(set((ROOT / "src").glob("*.F90")) |
                        set((ROOT / "utility").glob("*.F90")) |
                        set(HERE.rglob("*.F90")))
    modules = {}
    for path in candidates:
        for name in re.findall(r"^\s*module\s+(?!procedure\b)(\w+)",
                               path.read_text(), re.M | re.I):
            name = name.lower()
            if name in modules:
                raise ValueError(f"Duplicate module {name}")
            modules[name] = path
    ordered, visited, active = [], set(), set()

    def visit(path):
        if path in active:
            raise ValueError(f"Circular module dependency at {path}")
        if path in visited:
            return
        active.add(path)
        for name in re.findall(r"^\s*use\s+(\w+)", path.read_text(), re.M | re.I):
            name = name.lower()
            if name in {"machine", "netcdf", "mpi"}:
                continue
            if name not in modules:
                raise ValueError(f"Unresolved module {name} used by {path}")
            visit(modules[name])
        active.remove(path)
        visited.add(path)
        ordered.append(path)

    visit(HERE / "noahmp.F90")
    return ordered


if __name__ == "__main__":
    print(";".join(path.as_posix() for path in sources()))
