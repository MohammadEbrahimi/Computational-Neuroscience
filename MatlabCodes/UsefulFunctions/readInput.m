function []=readInput(fig_obj, eventDat,returnMap)

    key = eventDat.Key;
    switch key
        case 'leftarrow'
            newN = 2;
        case 'rightarrow'
            newN = 1;
        case 'uparrow'
            newN = 3;
        case 'downarrow'
            newN = 0;
    end

    returnMap('newN') = newN;
    close(fig_obj)
end