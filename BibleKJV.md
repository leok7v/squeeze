# "Bible Study"

Huffman + Deflate backreferences stats:

Compression Window: ```2^10 1024 bytes``` 

```4,436,173 -> 1,877,568 bytes 42.3% of "bible.txt"```

```nyt lit:80 len:28 pos:19 (count of NYT "Not Yet Transmitted" encoding points)```
Confirms that count of NYT encoding points is less or equal to the 
corresponding Huffman tree terminal nodes. 

```**bps:** 3.4 **H.lit:** 5.1 **H.pos: 3.5** **Bits** len: 4.3+0.1 pos: 3.8+6.4```

Average bits per symbol in the source tree is ```3.4```.

H.* is Shannon Entropy expressed in bits per symbol for
```lit``` huffman table that encodes literal bytes and first 30 base length indices.
```pos``` huffman table that encodes back reference positions indices.

**```Bits```** is average of Huffman bits + Extra bits

```lit: 19.5% of total 4,436,173 input bytes```
about 20% input bytes are passed thru to compressed stream
the rest it encoded as minimum length 3 backreferences.

```lit: 56.3% ref: 43.7% of total 1,539,161 encoding points```
The output compressed stream is composed of 56% of literal bytes
to 43% back references.


The length Huffman encoding table final composition:
```
len[ 3]: 37.5%  37.5% bits: 3
len[ 4]: 18.6%  56.0% bits: 4
len[ 5]: 13.0%  69.0% bits: 5
len[ 6]:  8.9%  77.9% bits: 5
len[ 7]:  5.8%  83.7% bits: 5
len[ 8]:  4.9%  88.6% bits: 6
len[ 9]:  3.7%  92.3% bits: 6
len[10]:  1.9%  94.2% bits: 6
len[11]:  1.2%  95.4% bits: 7
len[12]:  0.9%  96.3% bits: 7
len[13]:  0.7%  97.0% bits: 8
len[14]:  0.5%  97.5% bits: 8
len[15]:  0.4%  97.9% bits: 9
len[16]:  0.3%  98.3% bits: 9
len[17]:  0.3%  98.6% bits: 9
len[18]:  0.2%  98.8% bits: 9
len[19]:  0.2%  99.0% bits: 8
len[20]:  0.2%  99.1% bits: 8
len[21]:  0.1%  99.2% bits: 8
len[22]:  0.1%  99.3% bits: 8
```

The position Huffman encoding table final composition:
```
pos[ 1]:  0.0%   0.0% bits: 6
pos[ 2]:  0.1%   0.1% bits: 4
pos[ 3]:  0.1%   0.1% bits: 9
pos[ 4]:  0.2%   0.3% bits: 5
pos[ 5]:  0.2%   0.5% bits: 5
pos[ 6]:  0.7%   1.2% bits: 8
pos[ 7]:  1.0%   2.2% bits: 5
pos[ 8]:  2.1%   4.3% bits: 7
pos[ 9]:  2.2%   6.5% bits: 6
pos[10]:  4.2%  10.6% bits: 4
pos[11]:  3.5%  14.2% bits: 6
pos[12]:  7.0%  21.2% bits: 3
pos[13]:  5.9%  27.0% bits: 5
pos[14]: 10.0%  37.1% bits: 4
pos[15]:  8.1%  45.2% bits: 4
pos[16]: 12.9%  58.1% bits: 3
pos[17]: 10.5%  68.6% bits: 4
pos[18]: 17.2%  85.8% bits: 3
pos[19]: 14.2% 100.0% bits: 3
```

Full distance position frequency distribution is about this shape:
```
 #                                                              1.0
## #                                                           
## #                                                           
 # #                                                           
  ##                                                           
  ####                                                         
   ### #                                                       
     ###                                                        0.5
     ####                                                      
      #####                                                    
         ####                                                  
           #####                                               
            #######                                            
                ########## #                                   
                    #################                          
                          #  ########################### ## #  
                                           ####################
                                                               
                                                                0.0
     200                       512                         1023
```
see: 

and

commit:
