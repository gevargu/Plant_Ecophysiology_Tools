path1 = getDirectory("");// photos directory
list=getFileList(path1);// extract file names
path2 = getDirectory("");// blue channel directory
setBatchMode(true);// run in batch

for(i=0;i<list.length;i++){
    open(path1+list[i]);
    run("Rotate 90 Degrees Right");
    run("Rotate 90 Degrees Right");
    run("Split Channels");
    selectWindow(list[i]+" (green)");
    close();
    selectWindow(list[i]+" (red)");
    close();
    selectWindow(list[i]+" (blue)");
    saveAs("BMP",path2+list[i]+"_blue.bmp");
    close();
}

setBatchMode(false);
