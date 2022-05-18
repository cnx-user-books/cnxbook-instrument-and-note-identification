%basically takes the output of ProjectD and produces three beautiful graphs
%the first graph is the note rating over time of each note
%the second is a graph of the same data, but in a format similar to a musical score, but more similar to the scrolls used for player pianos
%the third is the "most likely" note over time, it is the easiest to read when only one note is being played at a time.
%the threshold argument is the threshold at which, if no note scores above, the sample is considered as silence and is graphed as -1 in the third graph.
function out = postProcessor(data,threshold) 
    figure; %draw note scores all in a big graph
    plot(data);
    title('note score over time');
    xlabel('window number');
    ylabel('score');
    legend('E0', 'F0', 'Fsharp0', 'G0', 'Gsharp0', 'A0', 'Asharp0', 'B0', 'C0', 'Csharp0', 'D0', 'Dsharp0', 'E1', 'F1', 'Fsharp1', 'G1', 'Gsharp1', 'A1', 'Asharp1', 'B1', 'C1', 'Csharp1', 'D1', 'Dsharp1', 'E2', 'F2', 'Fsharp2', 'G2', 'Gsharp2', 'A2', 'Asharp2', 'B2', 'C2', 'Csharp2', 'D2', 'Dsharp2');
    
    figure; %draw note scores in an imagesc
    set(gca, 'Ydir', 'normal');
    imagesc(data');
    title('note score over time');
    ylabel('window number');
    xlabel('note');
    
    figure; %draw most likely single note on a graph
    [A,B] = max(data');
    B = thresholdFunc(A,B,threshold);
    plot(B);
    title('best note score over time');
    xlabel('window number');
    ylabel('best note (1=A0)');

function out = thresholdFunc(A, B, threshold)
    for r = 1:length(A)
        if(A(r) < threshold)
           B(r) = -1;
        end
    end
    out = B;
    