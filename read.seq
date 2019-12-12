from random import choice


class Read:
    read: dict[int, int]
    count: int
    keys: list[int]
    size: int
    read_num: int
    mini_reads: dict[int, dict[int, int]]
    special_key: int
    rates: tuple[float, float]

    ## for normal reads, the read_dict doesnt include starts of reads
    def __init__(self: Read, read: dict[int, int], count: int, read_num: int):
        self.read = read
        self.count = count
        self.keys = sorted(self.read.keys())
        # self.genomic_region = None
        self.size = len(self.keys)
        self.read_num = read_num
        self.mini_reads = self.make_mini_reads()
        self.special_key = self.keys[0]
        self.rates = (0.5, 0.5)

    def __eq__(self: Read, other: Read):
        return self.read == other.read and self.read_num == other.read_num

    def __len__(self: Read):
        return len(self.read)

    def __str__(self: Read):
        return f'Read.{self.read}'

    def make_mini_reads(self: Read) -> dict[int, dict[int, int]]:
        mini_read_dict = dict[int, dict[int, int]]()
        mini_read = copy(self.read)
        for key in reversed(sorted(self.keys)):
            mini_read_dict[key] = mini_read
            mini_read = copy(mini_read)
            mini_read.pop(key)
        return mini_read_dict


def sample_from_reads(reads: list[Read]) -> list[Read]:
    new_reads = list[Read]()
    for read in reads:
        if read.size > 0:
            new_reads.append(read)
        else:
            choice = choice(read.keys)
            new_read = Read({choice: read.read[choice]}, 1, read.read_num)
            new_reads.append(new_read)
    return new_reads
