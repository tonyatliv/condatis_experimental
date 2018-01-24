import numpy as np
import landscape as ls

NORTH=1
SOUTH=2
EAST=3
WEST=4

# class AutoSourceTarget(ls.Landscape):
#     def __init__(self,size,edge,width):
#         rsx,rsy=size
#         rsx-=1
#         rsy-=1
#         if edge==NORTH:
#             x=(np.arange(rsx*width)%rsx).astype(np.int_)
#             y=((np.arange(rsx*width)/rsx)+rsy-width).astype(np.int_)
#         if edge==SOUTH:
#             x=(np.arange(rsx*width)%rsx).astype(np.int_)
#             y=(np.arange(rsx*width)/rsx).astype(np.int_)
#         if edge==EAST:
#             y=(np.arange(rsy*width)%rsy).astype(np.int_)
#             x=(np.arange(rsy*width)/rsy).astype(np.int_)
#         if edge==WEST:
#             y=(np.arange(rsy*width)%rsy).astype(np.int_)
#             x=((np.arange(rsy*width)/rsy)+rsx-width).astype(np.int_)
#         super(AutoSourceTarget, self).fromVecs(x,y)


def makeAutoSourceTarget(size,edge,width):
    rsx,rsy=size
    rsx-=1
    rsy-=1
    if edge==NORTH:
        x=(np.arange(rsx*width)%rsx).astype(np.int_)
        y=((np.arange(rsx*width)/rsx)+rsy-width).astype(np.int_)
    if edge==SOUTH:
        x=(np.arange(rsx*width)%rsx).astype(np.int_)
        y=(np.arange(rsx*width)/rsx).astype(np.int_)
    if edge==EAST:
        y=(np.arange(rsy*width)%rsy).astype(np.int_)
        x=(np.arange(rsy*width)/rsy).astype(np.int_)
    if edge==WEST:
        y=(np.arange(rsy*width)%rsy).astype(np.int_)
        x=((np.arange(rsy*width)/rsy)+rsx-width).astype(np.int_)
    return ls.Landscape.fromVecs(x,y)
    
