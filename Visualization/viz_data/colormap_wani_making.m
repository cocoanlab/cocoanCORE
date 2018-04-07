a = repmat(1:1000, 1000,1); imagesc(a);
colors = colormap(hsv);
colors(48:end,:) = [];
k=[ .65 .18 ]; col = k(1)*coljet+k(2); set(gcf,'colormap', flipud(col));

% There was some adjustment using colormapeditor after doing this.

% so here is my last colormap. 

col = [
    0.3843    0.1804    0.8314
    0.3216    0.1804    0.8314
    0.2627    0.1804    0.8314
    0.2000    0.1804    0.8314
    0.1804    0.2196    0.8314
    0.1804    0.2784    0.8314
    0.1804    0.3412    0.8314
    0.1804    0.4039    0.8314
    0.1804    0.4627    0.8314
    0.1804    0.5255    0.8314
    0.1804    0.5843    0.8314
    0.1804    0.6471    0.8314
    0.1804    0.6824    0.8314
    0.1804    0.7137    0.8314
    0.1804    0.7490    0.8314
    0.1804    0.7843    0.8314
    0.1804    0.8078    0.8314
    0.1804    0.8314    0.8314
    0.1804    0.8314    0.7059
    0.1804    0.8314    0.6471
    0.1804    0.8314    0.5843
    0.1804    0.8314    0.4627
    0.1804    0.8314    0.3412
    0.1804    0.8314    0.2196
    0.2627    0.8314    0.1804
    0.3843    0.8314    0.1804
    0.5059    0.8314    0.1804
    0.5725    0.8314    0.1804
    0.6392    0.8314    0.1804
    0.7059    0.8314    0.1804
    0.7725    0.8314    0.1804
    0.8157    0.8314    0.1804
    0.8314    0.8000    0.1804
    0.8314    0.7569    0.1804
    0.8314    0.7098    0.1804
    0.8314    0.6667    0.1804
    0.8314    0.6235    0.1804
    0.8314    0.5765    0.1804
    0.8314    0.5333    0.1804
    0.8314    0.4902    0.1804
    0.8314    0.4471    0.1804
    0.8314    0.4000    0.1804
    0.8314    0.3569    0.1804
    0.8314    0.3137    0.1804
    0.8314    0.2706    0.1804
    0.8314    0.2235    0.1804
    0.8314    0.1804    0.1804];

%% colors optimized for color-blind individuals (Wong, 2011, Nature Methods)

% colormap_blind_wani

cols_blind_8_names = {'1_Black', '2_Orange', '3_Sky blue', '4_Bluish green', '5_Yellow', '6_Blue', '7_Vermillion', '8_Reddish purple'}';
cols_blind_maps = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]/255;

save colormap_blind_wani cols_blind_8_names cols_blind_maps;


%% from http://colorbrewer2.org/index.php?type=diverging&scheme=RdYlBu&n=6

col6_from_colorbrewer = [
215 48 39
252 141 89
254 224 144
224 243 248
145 191 219
69 117 180] / 255;
col6_from_colorbrewer = flipud(col6_from_colorbrewer);

col11_from_colorbrewer = ...
    [158 1 66; 213 62 79; 244 109 67; 
    253 174 97; 254 224 139; 255 255 191; 
    230 245 152; 171 221 164; 102 194 165;
    50 136 189; 94 79 162]/255;
col11_from_colorbrewer = flipud(col11_from_colorbrewer);

save colormap_colorbrewer_wani col6_from_colorbrewer col11_from_colorbrewer;

%% ggplot colormap

col7_ggplot2 = [
    232 125 114
    190 154 51
    109 175 52
    86 188 151
    78 181 230
    159 143 248
    233 111 210]./255;

save colormap_ggplot2 col7_ggplot2;

col7_okabe_ito = [
    220 160 56
    109 179 228
    69 154 118
    239 226 98
    45 114 172
    199 100 38
    193 126 164]./255;
save solormap_okabe_ito col7_okabe_ito;

    
    