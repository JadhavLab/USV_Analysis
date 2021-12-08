classdef GoogleDrive < handle
    properties(Constant)
        URL_PREFIX="https://drive.google.com/file/d";
    end
    
    properties
        localFolder;
        pathKey;
    end
    
    methods
        function this=GoogleDrive(folderOnGoogleDrive)
            this.localFolder=fullfile(File.Home, 'Google Drive', folderOnGoogleDrive);
            if ~exist(this.localFolder, 'dir')
                msgError(['<html>Folder does not (yet) exist?' ...
                    Html.FileTree(this.localFolder) ...
                    '<br><br><center>???</center><hr></html>'],10)
            end
            this.pathKey=folderOnGoogleDrive;
        end
        
        function line=addFile(this, file)
            driveFile=fullfile(this.localFolder, file);
            if ~exist(driveFile, 'file')
                msgError(['<html>File does not (yet) exist?' ...
                    Html.FileTree(driveFile) '<hr></html>'],10)
                sz=0;
            else
                e=dir(driveFile);
                sz=e.bytes;
            end
            fileNoSpace=String.URLEncode(file, true);
            line=[fullfile(this.pathKey, fileNoSpace) '=' num2str(sz) ' '];
            googleLink=clipboard('paste');
            if startsWith(googleLink, GoogleDrive.URL_PREFIX)
            	line=[line googleLink];
            end
            if nargout==0
                clipboard('copy', line);
            end
        end
    end
    
    methods(Static)
        function HandleTooBig(dwl, modal)
            N=dwl.tooBigFiles.size;
            try
                if ishandle(gcf)
                    jw=Gui.JWindow(gcf);
                    jw.setAlwaysOnTop(false);
                end
            catch
            end
            warned=false;
            openedTarget=false;
            openedDownloads=false;
            if N<1
                return;
            end
            if nargin<2
                modal=true;
            end
            app=BasicMap.Global;
            fldrByFile=dwl.tooBigLocalFolderByRemoteFile;
            fldrByFldr=dwl.tooBigLocalFolderByRemoteFolder;
            nFldrs=dwl.tooBigLocalFolderByRemoteFile.size;
            pnl1=Gui.GridPanel([],nFldrs,1);
            it=fldrByFile.keySet.iterator;
            firstFldr='';
            while it.hasNext
                pnl1.add(BorderPnl(it.next));
            end
            ttl1='Google Drive automatic download limits..,.';
            if modal
                msgBox(Gui.Scroll(pnl1, 500, 550, app), ttl1);
            else
                msg(Gui.Scroll(pnl1, 500, 550, app), 10, 'center', ttl1 );
            end
            
            function bp=BorderPnl(fldr)
                if isempty(firstFldr)
                    firstFldr=fldr;
                end
                list=fldrByFile.get(fldr);
                bp=Gui.BorderPanel([],0,0);
                ttl=[num2str(list.size) ' file(s) too big for automatic download...'];
                Gui.SetTitledBorder(ttl, bp, 'bold', java.awt.Color.RED);
                p2=Gui.BorderPanel([],0,0);
                pnlFldrs=Gui.SetTitledBorder('Open folders?', p2, 'italic');
                p2.add(Gui.NewBtn(Html.WrapSmallBold('Local'),...
                    @(h,e)local(fldr, true)), 'North');
                p2.add(Gui.NewBtn(Html.WrapSmallBold('Downloads'),...
                    @(h,e)local(File.Downloads, false)), 'South');
                list2=fldrByFldr.get(fldr);
                if ~isempty(list2)
                    p2.add(Gui.NewBtn(Html.WrapSmallBold(...
                        ['Remote ' num2str(list2.size)]),...
                        @(h,e)remote(h, list2, true)), 'Center');
                end
                bp.add(Gui.FlowLeftPanel(4,0, pnlFldrs, ['<html>' ...
                    Html.FileTree(fldr, app,false) ...
                    '<hr><br></html>']), 'North');
                it=list.iterator;
                foci={'File name', 'Size', ['<html><b>' app.smallStart ...
                    'Browser<br>download?' app.smallEnd '</b></html>']};
                while it.hasNext
                    url=it.next;
                    [~, f, ext]=fileparts(char(url.driveKey));
                    foci{end+1}=String.URLDecode([f ext]);
                    delete(fullfile(fldr, foci{end})); % bad file
                    delete(fullfile(File.Downloads, foci{end})); % clear prior downloads
                    
                    foci{end+1}=String.encodeMb(url.size);
                    foci{end+1}=Gui.NewBtn(Html.WrapSmallBold(...
                        'Download' ), @(h,e)remote(h, ...
                        char(url.sharingUrl)),...
                        'Open browser to download file');
                end
                gbc=javaObjectEDT('java.awt.GridBagConstraints');
                anchors=[gbc.WEST gbc.EAST gbc.CENTER];
                center=Gui.GridBagPanel(0, 3, anchors, foci{:});
                bp.add(center, 'Center');
            end
            
            function remote(btn, url, isList)
                if ~warned
                end
                if nargin<3 || ~isList
                    web(url, '-browser');
                else
                    it=url.iterator;
                    while it.hasNext
                        web(char(it.next), '-browser');
                    end
                end
                if ~warned
                    jd=Gui.WindowAncestor(btn);
                    setAlwaysOnTopTimer(jd, 2, true, false);
                    MatBasics.RunLater(@(h,e)lazyOpen(jd), 2);
                end
            end
            
            function lazyOpen(jd)
                was=app.currentJavaWindow;
                app.currentJavaWindow=jd;
                app.currentJavaWindow.requestFocus;
                [~,f,ext]=fileparts(firstFldr);
                advice=msg(Html.WrapHr(['Once you have downloaded '...
                    '<br>file(s) drag and drop them<br>'...
                    'to the folder window titled<br>"<i>'...
                    f ext '</i>"']), 12, 'south+');
                app.currentJavaWindow=was;
                warned=true;
                if ~openedTarget
                    local(firstFldr, true);
                end
                if ~openedDownloads
                    local(File.Downloads, false);
                end
                setAlwaysOnTopTimer(advice, 2, true, false);
            end
            
            function local(fldr, isTarget)
                if isTarget
                    openedTarget=true;
                else
                    openedDownloads=true;
                end
                File.OpenFolderWindow(fullfile(fldr, '.'), '', false);
            end
        end
        
        function J=Test
            [~, ~, J]=WebDownload.Many({'omip044Labeled.csv',...
                'omip69Labeled.csv'});
        end
        
        function J=Test2
            disp('Expecting download window')
            [~, ~, J]=WebDownload.Many({'MC%203015036.fcs',...
                'MC%20303444.fcs'}, [], 'Samples/OMIP40Color');
        end
    end
end