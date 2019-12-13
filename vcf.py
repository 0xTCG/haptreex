# Modeled on https://raw.githubusercontent.com/jamescasbon/PyVCF/master/vcf/parser.py

# htslib integration
BAM_FPAIRED        =    1
BAM_FPROPER_PAIR   =    2
BAM_FUNMAP         =    4
BAM_FMUNMAP        =    8
BAM_FREVERSE       =   16
BAM_FMREVERSE      =   32
BAM_FREAD1         =   64
BAM_FREAD2         =  128
BAM_FSECONDARY     =  256
BAM_FQCFAIL        =  512
BAM_FDUP           = 1024
BAM_FSUPPLEMENTARY = 2048

from core.c_stubs import SAMHeaderTarget
from bio.align import CIGAR
from bio.locus import Locus

# This type must be consistent with htslib:
type _bam_core_t(tid: i32,
                 pos: i32,
                 bin: u16,
                 qual: u8,
                 l_extranul: u8,
                 flag: u16,
                 l_qname: u16,
                 n_cigar: u32,
                 l_qseq: i32,
                 mtid: i32,
                 mpos: i32,
                 isize: i32)

type SAMCore(_tid: i32, _pos: i32, _mapq: u8, _flag: u16, _mtid: i32, _mpos: i32, _isize: i32):
    def tid(self: SAMCore):
        return int(self._tid)

    def pos(self: SAMCore):
        return int(self._pos)

    def mapq(self: SAMCore):
        return int(self._mapq)

    def flag(self: SAMCore):
        return int(self._flag)

    def mtid(self: SAMCore):
        return int(self._mtid)

    def mpos(self: SAMCore):
        return int(self._mtid)

    def isize(self: SAMCore):
        return int(self._isize)

type SAMAux(s: ptr[u8]):
    def __bool__(self: SAMAux):
        return bool(self.s)

    @property
    def i(self: SAMAux):
        return _C.bam_aux2i(self.s)

    @property
    def f(self: SAMAux):
        return _C.bam_aux2f(self.s)

    @property
    def A(self: SAMAux):
        return _C.bam_aux2A(self.s)

    @property
    def B_len(self: SAMAux):
        return int(_C.bam_auxB_len(self.s))

    def B2i(self: SAMAux, idx: int):
        return _C.bam_auxB2i(self.s, u32(idx))

    def B2f(self: SAMAux, idx: int):
        return _C.bam_auxB2f(self.s, u32(idx))

class SAMRecord:
    _name: str 
    _read: seq 
    _qual: str 
    _cigar: CIGAR 
    _core: SAMCore 
    _aux: str
    
    def __init__(self: SAMRecord, aln: cobj) -> SAMRecord:
        name = _C.seq_hts_get_name(aln)
        read = _C.seq_hts_get_seq(aln)
        qual = _C.seq_hts_get_qual(aln)
        cigar = _C.seq_hts_get_cigar(aln)
        aux = _C.seq_hts_get_aux(aln)
        hts_core = ptr[_bam_core_t](aln)[0]
        core = SAMCore(hts_core.tid, hts_core.pos, hts_core.qual, hts_core.flag, hts_core.mtid, hts_core.mpos, hts_core.isize)
        return (name, read, qual, cigar, core, aux)

    @property
    def name(self: SAMRecord):
        return self._name

    @property
    def read(self: SAMRecord):
        return self._read

    @property
    def qual(self: SAMRecord):
        return self._qual

    @property
    def cigar(self: SAMRecord):
        return self._cigar

    @property
    def tid(self: SAMRecord):
        return self._core.tid()

    @property
    def pos(self: SAMRecord):
        return self._core.pos()

    @property
    def locus(self: SAMRecord):
        pos = self.pos
        return Locus(self.tid, -pos if self.reversed else pos)

    @property
    def mate_tid(self: SAMRecord):
        return self._core.mtid()

    @property
    def mate_pos(self: SAMRecord):
        return self._core.mpos()

    @property
    def mate_locus(self: SAMRecord):
        pos = self.mate_pos
        return Locus(self.mate_tid, -pos if self.mate_reversed else pos)

    @property
    def mapq(self: SAMRecord):
        return self._core.mapq()

    @property
    def insert_size(self: SAMRecord):
        return self._core.isize()

    @property
    def paired(self: SAMRecord):
        return self._core.flag() & BAM_FREVERSE != 0

    @property
    def proper_pair(self: SAMRecord):
        return self._core.flag() & BAM_FPROPER_PAIR != 0

    @property
    def unmapped(self: SAMRecord):
        return self._core.flag() & BAM_FUNMAP != 0

    @property
    def mate_unmapped(self: SAMRecord):
        return self._core.flag() & BAM_FMUNMAP != 0

    @property
    def reversed(self: SAMRecord):
        return self._core.flag() & BAM_FREVERSE != 0

    @property
    def mate_reversed(self: SAMRecord):
        return self._core.flag() & BAM_FMREVERSE != 0

    @property
    def read1(self: SAMRecord):
        return self._core.flag() & BAM_FREAD1 != 0

    @property
    def read2(self: SAMRecord):
        return self._core.flag() & BAM_FREAD2 != 0

    @property
    def secondary(self: SAMRecord):
        return self._core.flag() & BAM_FSECONDARY != 0

    @property
    def qc_fail(self: SAMRecord):
        return self._core.flag() & BAM_FQCFAIL != 0

    @property
    def duplicate(self: SAMRecord):
        return self._core.flag() & BAM_FDUP != 0

    @property
    def supplementary(self: SAMRecord):
        return self._core.flag() & BAM_FSUPPLEMENTARY != 0

    def aux(self: SAMRecord, tag: str):
        if len(tag) != 2:
            raise ValueError("SAM aux tags are two characters (got: " + tag + ")")
        return SAMAux(_C.seq_hts_aux_get(self._aux, tag))

extend SAMHeaderTarget:
    def __str__(self: SAMHeaderTarget):
        return self._name

    def __len__(self: SAMHeaderTarget):
        return self._len

class BAM:
    file: cobj
    idx: cobj
    hdr: cobj
    itr: cobj
    aln: cobj
    targets: list[SAMHeaderTarget]

    def _init(self: BAM, path: str, region: str):
        path_c_str, region_c_str = path.c_str(), region.c_str()

        file = _C.hts_open(path_c_str, "rb".c_str())
        if not file:
            raise IOError("file " + path + " could not be opened")

        idx = _C.sam_index_load(file, path_c_str)
        if not idx:
            _C.hts_close(file)
            raise IOError("unable to open BAM/CRAM index for " + path)

        hdr = _C.sam_hdr_read(file)
        itr = _C.sam_itr_querys(idx, hdr, region_c_str)

        if not hdr or not itr:
            _C.hts_close(file)
            raise IOError("unable to seek to region " + region + " in " + path)

        aln = _C.bam_init1()
        targets_array = _C.seq_hts_get_targets(hdr)
        targets = list[SAMHeaderTarget](targets_array, len(targets_array))

        self.file = file
        self.idx = idx
        self.hdr = hdr
        self.itr = itr
        self.aln = aln
        self.targets = targets

    def __init__(self: BAM, path: str):
        self._init(path, ".")

    def __init__(self: BAM, path: str, region: str):
        self._init(path, region)

    def _ensure_open(self: BAM):
        if not self.file:
            raise IOError("I/O operation on closed BAM/CRAM file")

    def _iter(self: BAM):
        self._ensure_open()
        while _C.seq_hts_sam_itr_next(self.file, self.itr, self.aln) >= 0:
            yield self.aln

        self.close()

    def __seqs__(self: BAM):
        for aln in self._iter():
            yield _C.seq_hts_get_seq(aln)

    def __iter__(self: BAM):
        for aln in self._iter():
            yield SAMRecord(aln)

    def close(self: BAM):
        if self.itr:
            _C.hts_itr_destroy(self.itr)

        if self.idx:
            _C.hts_idx_destroy(self.idx)

        if self.aln:
            _C.bam_destroy1(self.aln)

        if self.hdr:
            _C.bam_hdr_destroy(self.hdr)

        if self.file:
            _C.hts_close(self.file)

    def __enter__(self: BAM):
        pass

    def __exit__(self: BAM):
        self.close()

class SAM:
    file: cobj
    hdr: cobj
    aln: cobj
    targets: list[SAMHeaderTarget]

    def __init__(self: SAM, path: str):
        path_c_str = path.c_str()

        file = _C.hts_open(path_c_str, "r".c_str())
        if not file:
            raise IOError("file " + path + " could not be opened")

        hdr = _C.sam_hdr_read(file)
        aln = _C.bam_init1()
        targets_array = _C.seq_hts_get_targets(hdr)
        targets = list[SAMHeaderTarget](targets_array, len(targets_array))

        self.file = file
        self.hdr = hdr
        self.aln = aln
        self.targets = targets

    def _ensure_open(self: SAM):
        if not self.file:
            raise IOError("I/O operation on closed SAM file")

    def _iter(self: SAM):
        self._ensure_open()
        while True:
            status = int(_C.sam_read1(self.file, self.hdr, self.aln))
            if status >= 0:
                yield self.aln
            elif status == -1:
                break
            else:
                raise IOError("SAM read failed with status: " + str(status))

        self.close()

    def __seqs__(self: SAM):
        for aln in self._iter():
            yield _C.seq_hts_get_seq(aln)

    def __iter__(self: SAM):
        for aln in self._iter():
            yield SAMRecord(aln)

    def close(self: SAM):
        if self.aln:
            _C.bam_destroy1(self.aln)

        if self.hdr:
            _C.bam_hdr_destroy(self.hdr)

        if self.file:
            _C.hts_close(self.file)

        self.file = cobj()
        self.hdr = cobj()
        self.aln = cobj()

    def __enter__(self: SAM):
        pass

    def __exit__(self: SAM):
        self.close()

