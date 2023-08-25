# 修改

### hittable_list
- 實作kd-tree

### main.cc
- 讀入檔案物件後要先呼叫make_kd_tree()

### hittable.h
- 新增虛擬函數
```c++
  - get_kind(): 回傳物件類型
  - get_center(): 回傳物件重心
  - get_outer_bound(point3 *p): 回傳物件外圍點座標
```

### sphere.h
```c++
- 實作get_kind(): return 's'
- 實作get_center(): 回傳球心座標
- 實作get_outer_bound(point3 *p): 回傳球體座標中最大x、最小x、最大y、最小y、最大z、最小z的座標
```

### triangle.h
```c++
- 實作get_kind(): return 't'
- 實作get_center(): 回傳三角形重心座標
- 實作get_outer_bound(point3 *p): 回傳三角形三個點座標
```
<body>
