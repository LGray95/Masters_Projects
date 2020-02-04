def split_coords(coords):
    chr = coords.split(":")[0]
    start = int(coords.split(":")[1].split("-")[0])
    end = int(coords.split(":")[1].split("-")[1])
    return chr, start, end