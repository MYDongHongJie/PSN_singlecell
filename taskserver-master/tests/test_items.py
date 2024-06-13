import zipfile

import pytest

from taskserver.items import ProjectItem, Sfs2OBSItem


def test_list_tasks_dir(tmp_path):
    p = ProjectItem(root=str(tmp_path), code='project',
                    task_code='task', sub=True)
    p.workdir.mkdir(parents=True)
    p1 = p.workdir
    f1 = p1 / 'abc.txt'
    f1.write_text('h1')
    p11 = p1 / 'next'
    p11.mkdir()
    f11 = p11 / 'a.txt'
    f11.write_text('h2')

    detail = p.dir_stat()
    assert detail['size'] == 16  # 文件夹为4096, 两个文件各2,共4100
    assert detail['name'] == 'task'


class TestSfs2OBSItem:
    def test_to(self, mocker):
        obj = Sfs2OBSItem(source=['x'], op='haha')
        with pytest.raises(NotImplementedError):
            obj.to()
        func = mocker.patch('taskserver.items.Sfs2OBSItem.haha', create=True)
        obj.to()
        assert func.called

    def test_sfs2obs(self, mocker, tmp_path):
        p1 = tmp_path / 'p1'
        p1.mkdir()
        f1 = p1 / 'f1.txt'
        f1.write_text('f1')
        f2 = p1 / 'f2.txt'
        f2.write_text('f2')

        obj = Sfs2OBSItem(source=[str(p1)], target=[str(tmp_path)], op='sfs2obs')
        func1 = mocker.patch('taskserver.items.get_client',
                             **{'return_value.putFile.return_value': [('x', [('x1', {'reason': 'OK'})])]})
        obj.sfs2obs()
        assert func1.called
        assert not p1.exists()

    def test_obs2sfs(self, mocker, tmp_path):
        p1 = tmp_path / 'p1'
        p1.mkdir()
        p2 = tmp_path / 'p2'
        p2.mkdir()
        p2_z = p2 / 'p1.zip'
        with zipfile.ZipFile(p2_z, 'w', compression=zipfile.ZIP_DEFLATED) as handle:
            handle.writestr('f1.txt', 'f1')
            handle.writestr('f2.txt', 'f2')
            handle.writestr('d3/f1.txt', 'd3f1')
        obj = Sfs2OBSItem(source=[str(p1)], target=[str(p2)], op='obs2sfs')
        func1 = mocker.patch('taskserver.items.get_client',
                             **{'return_value.getObject.return_value.status': 200,
                                'return_value.deleteObject.return_value.status': 200})
        result = obj.obs2sfs()
        assert func1.called
        assert result['status'] == 0
        assert not p2_z.exists()
        assert len(list(p2.rglob('*'))) == 4
