"""
CNVision Gene Cache - Local Persistence Layer

Technologies: Python 3.10+, JSON (file format), file I/O

Core Functions:
- load_cache(): reads cached gene data from disk
- save_cache(): writes cache to persistent storage
- get_cached_gene(): retrieves gene from cache
- cache_gene(): adds gene to cache
"""

import json
import os

CACHE_FILE = "data/gene_cache.json"

def load_cache():
    """Load gene symbol cache from file."""
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            return json.load(f)
    return {}

def save_cache(cache):
    """Save gene symbol cache to file."""
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f)

def get_cached_gene(transcript_id, cache):
    """Get gene from cache."""
    return cache.get(transcript_id)

def cache_gene(transcript_id, gene_symbol, cache):
    """Add gene to cache."""
    cache[transcript_id] = gene_symbol