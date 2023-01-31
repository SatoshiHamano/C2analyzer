# coding: UTF-8

# モジュールのインポート
import os, tkinter, tkinter.filedialog, tkinter.messagebox, tkinter.ttk


def foo():
    print(t)
    print('Hello, %s!' % t.get())

# ファイル選択ダイアログの表示
root = tkinter.Tk()
root.withdraw()
fTyp = [("","*.py")]
iDir = os.path.abspath(os.path.dirname(__file__))
tkinter.messagebox.showinfo('○×プログラム','処理ファイルを選択してください！')
file = tkinter.filedialog.askopenfilename(filetypes = fTyp,initialdir = iDir)

# 処理ファイル名の出力
tkinter.messagebox.showinfo('○×プログラム',file)

root2 = tkinter.Tk()
root2.title("My First App")
frame1 = tkinter.ttk.Frame(root2)
label1 = tkinter.ttk.Label(frame1, text=file)

t = tkinter.StringVar()
entry1 = tkinter.ttk.Entry(frame1, textvariable=t)
button1 = tkinter.ttk.Button(frame1, text="OK", command=foo)

frame1.grid(row=0,column=0,sticky="nwse")
label1.grid(row=1,column=1,sticky="e")
entry1.grid(row=1,column=2,sticky="w")
button1.grid(row=2,column=2,sticky="w")

for child in frame1.winfo_children():
    child.grid_configure(padx=5, pady=5)

root2.mainloop()