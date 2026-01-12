import logging
from typing import Optional
from mpi4py import MPI

_comm: Optional[MPI.Comm] = None
_rank: int = -1
_size: int = -1
_initialized: bool = False

def init_mpi(comm: MPI.Comm = MPI.COMM_WORLD) -> None:
    """init MPI, only run onetime"""
    global _comm, _rank, _size, _initialized
    if _initialized:
        logging.warning("MPI has been initialized!!!")
        return
    
    _comm = comm
    _rank = comm.Get_rank()
    _size = comm.Get_size()
    _initialized = True
    logging.debug(f"MPI initialed: rank={_rank}, size={_size}")

def get_mpi_rank() -> int:
    """Get the local rank id"""
    if not _initialized:
        raise RuntimeError("Please run init_mpi() firstly!!!")
    return _rank

def get_mpi_size() -> int:
    """Get the total rank number"""
    if not _initialized:
        raise RuntimeError("Please run init_mpi() firstly!!!")
    return _size

def get_mpi_comm() -> MPI.Comm:
    """Get MPI communicator"""
    if not _initialized:
        raise RuntimeError("Please run init_mpi() firstly!!!")
    return _comm

def is_rank0() -> bool:
    """Check if current rank is rank 0"""
    return get_mpi_rank() == 0

def finalize_mpi() -> None:
    """finalize MPI"""
    global _comm, _rank, _size, _initialized
    if _initialized:
        _comm = None
        _rank = -1
        _size = -1
        _initialized = False