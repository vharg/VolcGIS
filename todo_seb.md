# TO-DO

- [ ] **Re-think the implementation of EPSG at eruption level**
    - I initially thought that using `EPSG:4326` would simplify things at a global scale, but it prevents computing geometry properties in meters. One option would be to switch to a Web Pseudo-Mercator `EPSG:3857` but I don't know how much distortion that would induce;
    - The other option is to use UTM Zones


- [ ] **Convert building footprints to points**
  - This seems required for saving to SpatiaLite