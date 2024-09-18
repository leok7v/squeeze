# Squeeze

zip like LZ77 + Huffman + Deflate compression

### Based on:

https://en.wikipedia.org/wiki/LZ77_and_LZ78

https://en.wikipedia.org/wiki/Adaptive_Huffman_coding

https://en.wikipedia.org/wiki/Deflate

### Goals:

* Simplicity (squeeze.h LoC: < 700).
* Ease of build and use (C99/C17/C23).
* Amalgamated into single header file library.
* No external dependencies.

### No goals:

* Performance (both CPU and memory).
* Existing archivers compatibility.
* Stream to stream encoding decoding.
* 16 and 32 bit CPU architectures.

### Overview:

The `squeeze` interface operates as a custom DEFLATE compression method, 
primarily relying on Huffman coding and LZ77 compression techniques. 
Here's a detailed breakdown of the core operations and data model:

### Data Model:
- **Huffman Trees**:
  - `squeeze` maintains two Huffman trees: one for literals/lengths (`lit`)
    and one for positions (`pos`). These trees are dynamically adjusted during
    compression and decompression.
  - Literal bytes (0-255) are represented as literal symbols, and the lengths
    (3-258) are assigned to symbols 257-285. A special "NYT" (Not Yet 
    Transmitted) symbol is used to handle newly encountered literals/positions.

- **LZ77 Backreferences**:
  - The compressor finds repeating patterns (length/distance pairs) in the 
    input stream and uses backreferences to encode them efficiently.
  - The distance for each backreference is encoded using the `pos` Huffman tree, 
    while the length is encoded using the `lit` Huffman tree.

### Key Methods:

#### 1. **Huffman Tree Operations**:
  - **`huffman_init`**: Initializes a Huffman tree with a predefined
    set of nodes. It ensures the tree starts in a balanced state where
    symbols are assigned equal frequencies.
  - **`huffman_insert`**: When a symbol (literal or length/position) is
    encountered for the first time, it is inserted into the Huffman tree,
    causing the tree to rebalance itself. This insertion involves potentially
    swapping siblings, adjusting tree paths, and recalculating symbol 
    frequencies.
  - **`huffman_inc_frequency`**: Updates the frequency of a symbol and 
    rebalances the Huffman tree accordingly. If a leaf node has a frequency 
    higher than its sibling, it may get promoted in the tree.

#### 2. **Compression Workflow**:
  - **`squeeze_compress`**: 
    - This method iterates over the input data. For each byte or sequence
      of bytes, it checks for backreferences to previous data within 
      a defined window (up to 32 KB).
    - If a backreference is found, the length and position are encoded 
      using the corresponding Huffman trees. If no backreference is found, 
      the byte is encoded as a literal.
    - After encoding, the Huffman trees are updated to reflect the increased 
      frequency of the symbol.
  
  - **`squeeze_write_huffman`**: Writes the Huffman-encoded value of a symbol
    to the output bitstream. After writing, it calls `huffman_inc_frequency` 
    to update the Huffman tree's structure.

#### 3. **Decompression Workflow**:
  - **`squeeze_decompress`**: Reverses the compression process. It reads 
    literals and backreferences from the bitstream, using the Huffman trees 
    to decode the values. The decoded data is then used to reconstruct 
    the original input stream.
  
  - **`squeeze_read_huffman`**: Reads a symbol from the bitstream, 
    following the Huffman tree structure. It returns the decoded symbol, 
    which can be either a literal byte or a length/position pair.

#### 4. **Bitstream Handling**:
  - **`squeeze_write_bit`/`squeeze_write_bits`**: These functions write
    individual bits or groups of bits to the bitstream. Compression
    algorithms often deal with partial byte data, so writing bits is
    essential to ensure the bitstream is packed efficiently.
  
  - **`squeeze_flush`**: Ensures that any pending bits in the bitstream
    are flushed to the output file or buffer.

### Huffman Encoding/Decoding:
- The `squeeze_len_base` and `squeeze_pos_base` arrays hold predefined
  base values for lengths and positions, respectively. These base values
  are combined with extra bits (stored in `squeeze_len_xb` and
  `squeeze_pos_xb`) to represent the full range of possible
  lengths and positions.
- When encoding a length or position, the compressor first writes the 
  corresponding Huffman code and then appends the extra bits required 
  to fully specify the value.
- Similarly, during decompression, the decoder reads the Huffman code 
  and any extra bits to reconstruct the full length or position value.

### Error Handling:
- The `error` field in the `squeeze_type` struct is used to track any 
  issues that arise during compression or decompression. 
  If an error occurs (e.g., out of memory, invalid input), the 
  compression/decompression process is halted.

### Theory of Operation Summary:
The `squeeze` interface provides an adaptive compression algorithm 
that dynamically adjusts its Huffman trees based on the input data. 
It uses LZ77 to find repeating patterns in the data and encodes 
them efficiently using backreferences. Huffman encoding is used to 
represent both literal bytes and length/position pairs compactly. 
By updating the Huffman trees as data is processed, the compressor 
adapts to the characteristics of the input data, ensuring that 
commonly occurring symbols are represented with fewer bits.

### Test materials:

Because Chinese texts are very compact comparing to e.g. the KJV bible
the Guttenberg License wording is stripped from the text files.

* See downloads.bat
