# Raytracer

### Usage
```
Usage: ./raytracer [OPTIONS]
1.
  --once                  render 1 frame
  --realtime              render continuously
2.
  --cornell-box           load test model
  --load     <file.obj>   load custom model
  --load-box <file.obj>   load custom model, inside Cornell Box
```

### Controls
```
 Move view:             Move light:             Rotate camera:

     W                      I                         ^
  A  S  D                J  K  L                   <  v  >

  Q --> up               U --> up
  E --> down             P --> down
```

# Rasteriser

### Usage
```
Usage: ./rasteriser [OPTIONS]
1.
  --once                  render 1 frame
  --realtime              render continuously
2.
  --cornell-box           load test model
  --load     <file.obj>   load custom model
  --load-box <file.obj>   load custom model, inside Cornell Box
```

### Controls
```
 Move view:             Move light:             Rotate camera:

     W                      I                       MOUSE  
  A  S  D                J  K  L                   

  Q --> up               U --> up
  E --> down             P --> down
```

Note: The `lucy.obj` model is not commited because its size is too big (2.45 GB)
