def get_bounds(bounds: int | tuple) -> bool:
    """
    The function outputs a range of values.
    
    Args:
        data: integer or tuple of integers.
    
    Return:
        two integers.
    """
    if type(bounds) == int:
        min_bound = 0
        max_bound = bounds
    
    else:
        min_bound = min(bounds)
        max_bound = max(bounds)

    return min_bound, max_bound


def filter_gc(sequence: str, gc_bounds: int | tuple=(0, 100)) -> bool:
    """
    The function filters readings by GC composition.
    
    Args:
        sequence: string, gc_bounds: integer or tuple of integers.
    
    Return:
        boolean result
    """
    gc = (sequence.count('G') + sequence.count('G')) / len(sequence)
    gc = gc * 100

    min_gc_bound, max_gc_bound = get_bounds(gc_bounds)
    
    return min_gc_bound <= gc <= max_gc_bound


def filter_length(sequence: str, length_bounds: int | tuple=(0, 2**32)) -> bool:
    """
    The function filters readings by length.
    
    Args:
        sequence: string, length_bounds: integer or tuple of integers.
    
    Return:
        boolean result
    """
    sequence_length = len(sequence)
    
    min_length_bound, max_length_bound = get_bounds(length_bounds)
    
    return  min_length_bound <= sequence_length <= max_length_bound


def filter_quality(sequence_quality: str, quality_threshold: int=0) -> bool:
    """
    The function filters readings by quality
    
    Args:
        sequence_quality: string, quality_threshold: integer or tuple of integers.
    
    Return:
        boolean result
    """
    
    sequence_quality_score = [(ord(chr_quality) - 33) for chr_quality in sequence_quality]
    sequence_quality_score = sum(sequence_quality_score) / len(sequence_quality_score)

    return quality_threshold <= sequence_quality_score