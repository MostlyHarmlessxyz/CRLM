"""Read only HDF5 metadata/obs from the public atlas using bounded HTTP ranges.

Never reads the expression X matrix. Requires 206 Content-Range responses and
refuses whole-object downloads. At most 64 MiB is fetched in any invocation.
"""
import io
import json
import urllib.request
from pathlib import Path
import h5py

ROOT = Path(__file__).resolve().parent
URL = "https://cellgeni.cog.sanger.ac.uk/pan-immune/global.h5ad"
SIZE = 5131003648
CACHE = ROOT / "geo/normal_atlas_range_cache"
CACHE.mkdir(exist_ok=True)

class RangeReader(io.RawIOBase):
    block = 1048576
    def __init__(self):
        self.pos = 0
        self.transferred = 0
        self.requests = 0
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else self.pos + offset if whence == 1 else SIZE + offset
        if self.pos < 0: raise ValueError("Negative seek")
        return self.pos
    def read(self, n=-1):
        if n < 0: raise ValueError("Unbounded reads prohibited")
        n = min(n, SIZE - self.pos)
        chunks = []
        while n > 0:
            index = self.pos // self.block
            target = CACHE / f"block_{index}.bin"
            if not target.exists():
                start = index * self.block
                stop = min(start + self.block, SIZE) - 1
                if self.transferred + (stop-start+1) > 64 * 1048576:
                    raise ValueError("Metadata transfer cap reached")
                request = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{stop}"})
                for attempt in range(3):
                    try:
                        with urllib.request.urlopen(request, timeout=40) as response:
                            if response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{stop}/{SIZE}":
                                raise ValueError("Server ignored exact range; refusing whole file")
                            body = response.read(stop-start+2)
                        break
                    except OSError:
                        if attempt == 2: raise
                if len(body) != stop-start+1: raise ValueError("Wrong range response size")
                target.write_bytes(body)
                self.transferred += len(body)
                self.requests += 1
                if self.requests % 5 == 0: print("metadata ranges", self.requests, "bytes", self.transferred, flush=True)
            body = target.read_bytes()
            offset = self.pos % self.block
            take = min(n, len(body)-offset)
            if take <= 0: raise ValueError("Empty cached range")
            chunks.append(body[offset:offset+take])
            self.pos += take
            n -= take
        return b"".join(chunks)
    def readinto(self, b):
        data = self.read(len(b)); b[:len(data)] = data; return len(data)

def decoded(values):
    return [x.decode() if isinstance(x,bytes) else x.item() if hasattr(x,'item') else x for x in values]

reader = RangeReader()
with h5py.File(reader, "r") as atlas:
    print("ROOT", list(atlas.keys()), flush=True)
    print("OBS", list(atlas["obs"].keys()), flush=True)
    for key, value in atlas["obs"].items():
        print(key, "group" if isinstance(value,h5py.Group) else str(value.shape)+" "+str(value.dtype), flush=True)
    selected = {}
    for key in atlas["obs"].keys():
        if any(term in key.lower() for term in ["donor", "tissue", "organ", "sample", "annotation", "celltypist", "cell_type", "celltype", "cell_category"]):
            obj = atlas["obs"][key]
            if isinstance(obj,h5py.Group) and "codes" in obj and "categories" in obj:
                codes = obj["codes"][:]
                categories = decoded(obj["categories"][:])
                selected[key] = [categories[i] if i >= 0 else "unknown" for i in codes]
            elif isinstance(obj,h5py.Dataset):
                values = obj[:]
                if "__categories" in atlas["obs"] and key in atlas["obs/__categories"]:
                    categories = decoded(atlas["obs/__categories"][key][:])
                    selected[key] = [categories[i] if i >= 0 else "unknown" for i in values]
                else: selected[key] = decoded(values)
    (ROOT / "geo/E-MTAB-11536_obs_selected.json").write_text(json.dumps(selected))
    print("selected_columns", list(selected), "rows", {k:len(v)for k,v in selected.items()}, flush=True)
print("network_bytes", reader.transferred, "requests", reader.requests, flush=True)
